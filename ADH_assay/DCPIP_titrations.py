from opentrons import protocol_api
from opentrons.protocol_api import ALL, PARTIAL_COLUMN, SINGLE
from opentrons.types import Point
import time
import sys
import math
import random
import subprocess

metadata = {
    'protocolName': 'DCPIP Standard Curve for Extinction Coefficient',
    'author': 'Sofia Terenziani',
    'description': 'This protocol performs a 40/60 DCPIP serial dilution to create a standard curve and determine the extinction coefficient of DCPIP. The plate will be run at 3 different pHs (6,7,8)'
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
    global tips_1000, plate, buffer, pipette, trash, dcpip
    tips_1000 = protocol.load_labware('opentrons_flex_96_tiprack_1000ul', 'B1')
    plate = protocol.load_labware('corning_384_wellplate_112ul_flat', 'C2')
    buffer= protocol.load_labware('nest_12_reservoir_15ml', 'C3')
    pipette = protocol.load_instrument('flex_96channel_1000', 'right', tip_racks=[tips_1000])
    dcpip = protocol.load_labware('greiner_96_wellplate_323ul', 'B3')
    trash = protocol.load_trash_bin ('D1')

def define_liquids(protocol):
    global dcpip, buffer_wells, dcpip_liquid
    buffer_wells = [buffer.wells_by_name()[f'A{i}'] for i in range(1, 13)] # Wells 1,2,3,4 pH6 - wells 5,6,7,8 pH7 - wells 9,10,11,12 pH8

    buffer_liquid_light = protocol.define_liquid(
        name="Buffer pH6",
        description="150 mM NaCl, 100 mM HEPES, pH 6",
        display_color="#ADD8E6")
    buffer_liquid_medium = protocol.define_liquid(
        name="Buffer pH7",
        description="150 mM NaCl, 100 mM HEPES, pH 7",
        display_color="#6495ED")
    buffer_liquid_dark = protocol.define_liquid(
        name="Buffer pH8",
        description="150 mM NaCl, 100 mM HEPES, pH 8",
        display_color="#00008B")
    
    dcpip = dcpip.rows()[0] # 500 uM initial concentration of DCPIP (150uM concentration in 384 well plate rxns)
    dcpip_liquid = protocol.define_liquid(
        name="DCPIP",
        description="500 µM DCPIP solution",
        display_color="#006400")

    for well in buffer_wells[0:4]:
        well.load_liquid(liquid=buffer_liquid_light, volume=15000)
 
    for well in buffer_wells[4:8]:
        well.load_liquid(liquid=buffer_liquid_medium, volume=15000)
 
    for well in buffer_wells[8:12]:
        well.load_liquid(liquid=buffer_liquid_dark, volume=15000)

    for well in dcpip:
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
    for col_idx in range(2):
        pickup_tips('row', protocol)
        col = plate.columns()[col_idx]
        buffer_source = buffer_wells[col_idx % len(buffer_wells)]
        col_wells = col[1:16]
        pipette.distribute(60, buffer_source, col_wells, new_tip='never', mix_before=(3,30), disposal_volume=10)
        pipette.drop_tip() 
 
def add_and_titrate_dcpip(protocol):
    pickup_tips('row', protocol)
    pipette.aspirate(110,dcpip)
    pipette.dispense(100, plate.rows()[0])
    pipette.drop_tip()

    for col_idx in range(2):
        pickup_tips('row', protocol)
        col = plate.columns()[col_idx]
        col_wells = col[0:15]
        pipette.distribute(40, dcpip, col_wells, new_tip='never', mix=(3,40), disposal_volume=40)
        pipette.aspirate(20, plate.rows()[15][0])
        pipette.drop_tip()

    pickup_tips('row', protocol)