from opentrons import protocol_api
from opentrons.protocol_api import COLUMN, ROW, SINGLE, ALL
from opentrons.types import Point

metadata = {
    'protocolName': 'PQQ-ADH Assay',
    'author': 'Sofia Terenziani, Shawn Laursen',
    'description': '''
    Adds buffer, PQQ, DCPIP and PMS to each well of 384 plate
    Adds each of the 23 metals mixtures + buffer control
    Adds each of the 15 alcohol mixtures + buffer control
    Adds PQQ-ADH mixtures to each well of 384 plate'''}

requirements = {'robotType': 'Flex','apiLevel': '2.29'}

# Runtime options

def add_parameters(parameters):
    parameters.add_bool(
        variable_name="USE_TEMP_MODULE",
        display_name="Temp module in slot C1?",
        description="Temp module in slot C1?",
        default=False)
    
def run(protocol):
    protocol.set_rail_lights(True)
    setup(protocol)
    define_liquids(protocol)
    add_buffer(protocol)
    add_metals(protocol)
    add_alcohols(protocol)
    add_controls(protocol)
    add_PQQ_ADH(protocol)
    protocol.set_rail_lights(False)

def setup(protocol):
    # Load modules and labware
    global tips_96_buffer, tips_96_enzyme, tips_rows, tips_columns, plate, metals, alcohols, enzyme, trash, buff_pqq_dcpip_pms, pipette, temp_module
    tips_96_buffer = protocol.load_labware('opentrons_flex_96_tiprack_1000ul', 'B3', adapter='opentrons_flex_96_tiprack_adapter')
    tips_96_enzyme = protocol.load_labware('opentrons_flex_96_tiprack_1000ul', 'A3', adapter='opentrons_flex_96_tiprack_adapter')
    tips_rows = protocol.load_labware('opentrons_flex_96_tiprack_1000ul', 'B1')
    tips_columns = protocol.load_labware('opentrons_flex_96_tiprack_1000ul', 'A2')
    # Optionally load the temperature module in C1 or use plain labware in C1
    if protocol.params.USE_TEMP_MODULE is True:
        temp_module = protocol.load_module('temperature module gen2', 'C1')
    else:
        temp_module = None

    # Labware
    plate = protocol.load_labware('corning_384_wellplate_112ul_flat', 'C2')
    metals = protocol.load_labware('greiner_96_wellplate_323ul', 'D2')
    alcohols = protocol.load_labware('greiner_96_wellplate_323ul', 'C3')
    enzyme = protocol.load_labware('nest_1_reservoir_195ml', 'B2')
    trash = protocol.load_trash_bin ('D1')
    # Load the buffer/PQQ/DCPIP/PMS reservoir either on the temp module or directly in C1
    if temp_module is not None:
        buff_pqq_dcpip_pms = temp_module.load_labware('nest_1_reservoir_195ml')
    else:
        buff_pqq_dcpip_pms = protocol.load_labware('nest_1_reservoir_195ml', 'C1')
    pipette = protocol.load_instrument('flex_96channel_1000')

    #volumes
    global buffer_volume, metals_volume, alcohols_volume, enzyme_volume
    rxn_vol = 60
    buffer_volume = rxn_vol/4
    metals_volume = rxn_vol/8
    alcohols_volume = rxn_vol/8
    enzyme_volume = rxn_vol/2

def define_liquids(protocol):
    buffer_liquid = protocol.define_liquid(
        name="Buffer/PQQ/DCPIP/PMS Mix",
        description="Buffer mixture containing PQQ, DCPIP, and PMS",
        display_color="#50C878")
    
    metals_liquid = protocol.define_liquid(
        name="Metal Mixtures",
        description="23 different metal mixtures",
        display_color="#FFD700")
    
    alcohols_liquid = protocol.define_liquid(
        name="Alcohol Mixtures",
        description="15 different alcohol mixtures",
        display_color="#FF6B6B")
    
    water_liquid = protocol.define_liquid(
        name="Water Control",
        description="Water control for negative control wells",
        display_color="#00BFFF")
    
    tcep_liquid = protocol.define_liquid(
        name="TCEP Control",
        description="TCEP control for stability checks",
        display_color="#8A2BE2")
    
    enzyme_liquid = protocol.define_liquid(
        name="PQQ-ADH Enzyme",
        description="15 PQQ-ADH enzyme mixtures",
        display_color="#4169E1")
    
    buff_pqq_dcpip_pms['A1'].load_liquid(liquid=buffer_liquid,volume=195000)
    
    for i in range(2):
        for j in range(12): metals.rows()[i][j].load_liquid(liquid=metals_liquid,volume=300)
            
    for i in range(16):
        alcohols.wells()[i].load_liquid(liquid=alcohols_liquid,volume=300)

    alcohols.rows()[7][1].load_liquid(liquid=water_liquid, volume=300)
    alcohols.rows()[0][2].load_liquid(liquid=tcep_liquid, volume=300)
    
    enzyme['A1'].load_liquid(liquid=enzyme_liquid,volume=195000)

def pickup_tips(layout, protocol):
    if layout == 'column':
        pipette.configure_nozzle_layout(style=protocol_api.COLUMN,start="A12", tip_racks=[tips_columns])   
    elif layout == 'row':
        pipette.configure_nozzle_layout(style=protocol_api.ROW,start="H1",tip_racks=[tips_rows])
    elif layout == 'single':
        pipette.configure_nozzle_layout(style=protocol_api.SINGLE,start="A1",tip_racks=[tips_rows])
    elif layout == 'all':
        pipette.configure_nozzle_layout(style=protocol_api.ALL,start="A1",tip_racks=[tips_96_buffer, tips_96_enzyme])
    pipette.pick_up_tip()

def add_buffer(protocol):
    pickup_tips('all', protocol)
    destinations = [plate.wells_by_name()[well] for well in ['A1', 'A2', 'B2', 'B1']]
    source = buff_pqq_dcpip_pms.wells_by_name()['A1']
    total_volume = buffer_volume * len(destinations) + 20
    pipette.mix(2, 80, source)
    pipette.aspirate(total_volume, source.bottom(1))
    for dest in destinations:
        pipette.dispense(buffer_volume, dest.bottom(1))
        pipette.touch_tip(location=dest, v_offset=-4, speed=20)
    pipette.blow_out(source.top())
    pipette.return_tip()

def add_metals(protocol):
    for i in range(2):
        pickup_tips('row', protocol)
        destinations = [plate.rows()[row][i] for row in range(15)]
        source = metals.rows()[i][0]
        total_volume = metals_volume * len(destinations) + 20
        pipette.mix(2, 100, source)
        pipette.aspirate(total_volume, source.bottom(1))
        for dest in destinations:
            pipette.dispense(metals_volume, dest.top())
            pipette.touch_tip(location=dest, v_offset=-3, speed=20)
        pipette.blow_out(source.top())
        pipette.drop_tip()

def add_alcohols(protocol):
    for j in range(2):
        pickup_tips('column', protocol)
        destinations = [plate.rows()[j][column] for column in range(24)]
        source = alcohols.rows()[0][j]
        total_volume = (alcohols_volume * len(destinations)) + 20
        pipette.mix(2, 100, source)
        pipette.aspirate(total_volume, source.bottom(1))
        for dest in destinations:
            pipette.dispense(alcohols_volume, dest.top())
            pipette.touch_tip(location=dest, v_offset=-2, speed=20)
        pipette.drop_tip()

def add_controls(protocol):
    pickup_tips('single', protocol)
    water = alcohols.rows()[7][2]
    tcep = alcohols.rows()[0][2]

    water_dest = [plate.rows()[15][column] for column in range(12)]
    tcep_dest = [plate.rows()[15][column] for column in range(12, 24)]

    for source, destinations in [(water, water_dest), (tcep, tcep_dest)]:
        total_volume = (alcohols_volume * len(destinations)) + 20
        pipette.mix(2, 100, source)
        pipette.aspirate(total_volume, source.bottom(1))
        for dest in destinations:
            pipette.dispense(alcohols_volume, dest.top())
            pipette.touch_tip(location=dest, v_offset=-2, speed=20)
            pipette.blow_out(source.top())
    pipette.drop_tip()


def add_PQQ_ADH(protocol):
    pickup_tips('all', protocol)
    destination_wells = [plate.wells_by_name()[well] for well in ['A1', 'A2', 'B2', 'B1']]
    source = enzyme['A1']
    total_volume = enzyme_volume * len(destination_wells) + 20
    pipette.aspirate(total_volume+20, source.bottom(1))
    pipette.dispense(10, source)
    for dest in destination_wells:
        pipette.dispense(enzyme_volume, dest.top(z=-3))
        pipette.touch_tip(dest, v_offset=-1, speed=20)
    pipette.blow_out(source.top())
    pipette.return_tip()