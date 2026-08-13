from opentrons import protocol_api
from opentrons.protocol_api import COLUMN, ROW, SINGLE, ALL
from opentrons.types import Point

metadata = {
    'protocolName': 'DCPIP Standard Curve for Extinction Coefficient',
    'author': 'Sofia Terenziani',
    'description': 'This protocol performs a 1:6 DCPIP serial dilution to create a standard curve and determine the extinction coefficient of DCPIP.'
}

requirements = {'robotType': 'Flex','apiLevel': '2.29'}

def run (protocol):
    protocol.set_rail_lights(True)
    setup(protocol)
    define_liquids(protocol)
    add_buffer(protocol)
    add_and_titrate_dcpip(protocol)
    protocol.set_rail_lights(False)

def setup(protocol):
    # Load labware
    global tips_1000, plate, reservoir, pipette, trash
    tips_1000 = protocol.load_labware('opentrons_flex_96_tiprack_1000ul', 'B2')
    plate = protocol.load_labware('corning_384_wellplate_112ul_flat', 'C2')
    reservoir= protocol.load_labware('nest_12_reservoir_15ml', 'C3')
    pipette = protocol.load_instrument('flex_96channel_1000', 'right', tip_racks=[tips_1000])
    trash = protocol.load_trash_bin ('D1')

def define_liquids(protocol):
    global dcpip, buffer_wells, buffer_liquid, dcpip_liquid
    buffer_wells = [reservoir.wells_by_name()[well] for well in ['A1', 'A2']] # 150 mM NaCl, 100 mM HEPES
    buffer_liquid = protocol.define_liquid(name="Buffer",description="150 mM NaCl, 100 mM HEPES",display_color="#ADD8E6")
    dcpip = reservoir.wells_by_name()['A3'] # 500 uM initial concentration of DCPIP (150uM final concentration in 384 well plate)
    dcpip_liquid = protocol.define_liquid(name="DCPIP",description="500 µM DCPIP solution",display_color="#006400")

    for well in buffer_wells:
        well.load_liquid(liquid=buffer_liquid, volume=15000)
    for well in [dcpip]:
        well.load_liquid(liquid=dcpip_liquid, volume=10000)

def pickup_tips(layout, protocol):
    if layout == 'column':
        pipette.configure_nozzle_layout(style=protocol_api.COLUMN,start="A12", tip_racks=[tips_1000])
    elif layout == 'row':
            pipette.configure_nozzle_layout(style=protocol_api.ROW,start="H1",tip_racks=[tips_1000])
    elif layout == 'single':
            pipette.configure_nozzle_layout(style=protocol_api.SINGLE,start="A1",tip_racks=[tips_1000])
    elif layout == 'all':
            pipette.configure_nozzle_layout(style=protocol_api.ALL,start="A1",tip_racks=[tips_1000])
    pipette.pick_up_tip()

def add_buffer(protocol):
    global add_buffer,buffer_wells,buffer_liquid
    rxn_vol = 60
    pickup_tips('column', protocol)
    row_a_wells = plate.rows()[0][1:24]
    row_b_wells = plate.rows()[1][1:24]
    destination_wells = row_a_wells + row_b_wells
    pipette.distribute(rxn_vol, buffer_wells[0], row_a_wells,new_tip='never', mix_before=(3, rxn_vol/2),disposal_volume=10)
    pipette.distribute(rxn_vol, buffer_wells[1], row_b_wells,new_tip='never', mix_before=(3, rxn_vol/2),disposal_volume=10)
    pipette.drop_tip()


def add_and_titrate_dcpip(protocol):
    global dcpip
    rxn_vol = 60
    dilution_factor = 6
    num_dilutions = 23
    transfer_vol = rxn_vol / dilution_factor
    row_a_wells = plate.rows()[0][0:24] 
    row_b_wells = plate.rows()[1][0:24]
    
    pickup_tips('column', protocol)
    pipette.transfer(rxn_vol + transfer_vol,dcpip, row_a_wells[0], mix_before=(3, rxn_vol/2), new_tip='never',mix_after=(3, rxn_vol/2))
    pipette.transfer(transfer_vol, row_a_wells[0:num_dilutions-1], row_a_wells[1:num_dilutions], new_tip='never', mix_after=(3, rxn_vol/2))
    pipette.transfer(rxn_vol + transfer_vol,dcpip, row_b_wells[0],mix_before=(3, rxn_vol/2),new_tip='never',mix_after=(3, rxn_vol/2))
    pipette.transfer(transfer_vol,row_b_wells[0:num_dilutions-1],row_b_wells[1:num_dilutions],new_tip='never',mix_after=(3, rxn_vol/2))
    pipette.drop_tip()