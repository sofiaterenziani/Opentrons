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

def run(protocol):
    protocol.set_rail_lights(True)
    setup(protocol)
    define_liquids(protocol)
    add_buffer(protocol)
    add_metals(protocol)
    add_alcohols(protocol)
    add_PQQ_ADH(protocol)
    protocol.set_rail_lights(False)

def setup(protocol):
    # Load modules and labware
    global tips_96, tips_rows, tips_columns, plate, metals, alcohols, enzyme, trash, buff_pqq_dcpip_pms, pipette
    tips_96 = protocol.load_labware('opentrons_flex_96_tiprack_200ul', 'B1', adapter='opentrons_flex_96_tiprack_adapter')
    tips_rows = protocol.load_labware('opentrons_flex_96_tiprack_200ul', 'B3')
    tips_columns = protocol.load_labware('opentrons_flex_96_tiprack_200ul', 'D3')

    # Labware
    plate = protocol.load_labware('corning_384_wellplate_112ul_flat', 'C2')
    metals = protocol.load_labware('greiner_96_wellplate_323ul', 'D2')
    alcohols = protocol.load_labware('greiner_96_wellplate_323ul', 'C3')
    enzyme = protocol.load_labware('nest_1_reservoir_195ml', 'B2')
    trash = protocol.load_trash_bin ('D1')
    buff_pqq_dcpip_pms = protocol.load_labware('nest_1_reservoir_195ml', 'C1')
    pipette = protocol.load_instrument('flex_96channel_1000')

    #volumes
    global buffer_volume, metals_volume, alcohols_volume, enzyme_volume
    rxn_vol = 40
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
    
    enzyme_liquid = protocol.define_liquid(
        name="PQQ-ADH Enzyme",
        description="15 PQQ-ADH enzyme mixtures",
        display_color="#4169E1")
    
    buff_pqq_dcpip_pms['A1'].load_liquid(liquid=buffer_liquid,volume=195000)
    
    for i in range(2):
        for j in range(12): metals.rows()[i][j].load_liquid(liquid=metals_liquid,volume=300)
            
    for i in range(16):
        alcohols.wells()[i].load_liquid(liquid=alcohols_liquid,volume=300)
    
    enzyme['A1'].load_liquid(liquid=enzyme_liquid,volume=195000)

def pickup_tips(layout, protocol):
    if layout == 'column':
        pipette.configure_nozzle_layout(style=protocol_api.COLUMN,start="A12", tip_racks=[tips_columns])   
    elif layout == 'row':
        pipette.configure_nozzle_layout(style=protocol_api.ROW,start="H1",tip_racks=[tips_rows])
    elif layout == 'all':
        pipette.configure_nozzle_layout(style=protocol_api.ALL,start="A1",tip_racks=[tips_96])
    pipette.pick_up_tip()

def add_buffer(protocol):
    pickup_tips('all', protocol)
    destinations = [plate.wells_by_name()[well] for well in ['A1', 'A2', 'B2', 'B1']]
    pipette.distribute(buffer_volume, buff_pqq_dcpip_pms.wells_by_name()['A1'], destinations, new_tip='never', disposal_volume=20)
    pipette.return_tip()

def add_metals(protocol):
    for i in range(2):
        pickup_tips('row', protocol)
        destinations = [plate.rows()[row][i].top().move(Point(x=1, y=0, z=-1)) for row in range(16)]
        pipette.distribute(metals_volume, metals.rows()[i][0], destinations, new_tip='never', disposal_volume=20)
        pipette.drop_tip()

def add_alcohols(protocol):
    for j in range(2):
        pickup_tips('column', protocol)
        destinations = [plate.rows()[j][column].top().move(Point(x=-1, y=0, z=-1)) for column in range(24)]
        pipette.distribute(alcohols_volume, alcohols.rows()[0][j], destinations, new_tip='never', disposal_volume=20)
        pipette.drop_tip()

def add_PQQ_ADH(protocol):
    tips_96.reset()
    pickup_tips('all', protocol)
    destinations = [plate.wells_by_name()[well].top().move(Point(x=0, y=1, z=-1)) for well in ['A1', 'A2', 'B2', 'B1']]
    pipette.distribute(enzyme_volume, enzyme['A1'], destinations, new_tip='never', disposal_volume=20)
    pipette.return_tip()