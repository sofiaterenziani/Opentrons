from opentrons import protocol_api
from opentrons.protocol_api import ALL, PARTIAL_COLUMN, SINGLE
from opentrons.types import Point
import time
import sys
import math
import random
import subprocess

metadata = {
    'protocolName': 'DCPIP titrations with OT2',
    'author': 'Sofia Terenziani',
    'description': 'Protocol for DCPIP titrations using the Opentrons OT-2 robot, starting at 500uM DCPIP, and testing at pH 6,7,8',
}

requirements = {
    'robotType': 'OT-2', 'apiLevel': '2.16'
}

def run (protocol):
    protocol.set_rail_lights(True)
    setup(protocol)
    define_liquids(protocol)
    add_20uL_buffer(protocol)
    add_and_titrate_dcpip(protocol)
    #add_20uL_buffer_and_mix(protocol)
    protocol.set_rail_lights(False)

def setup(protocol):
    # Load Labware
    global tips_300, plate, buffer, p300m, dcpip
    tips_300 = protocol.load_labware('opentrons_96_tiprack_300ul', '1')
    plate = protocol.load_labware('corning_384_wellplate_112ul_flat', '2')
    buffer = protocol.load_labware('nest_12_reservoir_15ml', '3')
    p300m = protocol.load_instrument('p300_multi_gen2', 'left', tip_racks=[tips_300])
    dcpip = protocol.load_labware('greiner_96_wellplate_323ul', '4')

    # Reagents
    global buffer_pH6, buffer_pH7, buffer_pH8, dcpip_pH6, dcpip_pH7, dcpip_pH8, liquid_waste
    buffer_pH6 = [buffer.wells_by_name()[f'A{i}'] for i in range(1, 5)]
    buffer_pH7 = [buffer.wells_by_name()[f'A{i}'] for i in range(5, 9)]
    buffer_pH8 = [buffer.wells_by_name()[f'A{i}'] for i in range(9, 13)]
    dcpip_pH6 = dcpip.columns()[0]
    dcpip_pH7 = dcpip.columns()[1]
    dcpip_pH8 = dcpip.columns()[2]
    liquid_waste = dcpip.columns()[11][0]

    # Define liquids

def define_liquids(protocol):
    """Add color-coded reagents and their starting volumes to the visualization."""
    buffer_pH6_liquid = protocol.define_liquid(
        name='Buffer pH 6', description='MES Buffer at pH 6', display_color='#ADD8E6')
    buffer_pH7_liquid = protocol.define_liquid(
        name='Buffer pH 7', description='HEPES Buffer at pH 7', display_color='#5B9BD5')
    buffer_pH8_liquid = protocol.define_liquid(
        name='Buffer pH 8', description='HEPES Buffer at pH 8', display_color='#0B3D91')
    dcpip_pH6_liquid = protocol.define_liquid(
        name='DCPIP pH 6', description='DCPIP titration solution at pH 6', display_color="#00A650")
    dcpip_pH7_liquid = protocol.define_liquid(
        name='DCPIP pH 7', description='DCPIP titration solution at pH 7', display_color="#00A6509E")
    dcpip_pH8_liquid = protocol.define_liquid(
        name='DCPIP pH 8', description='DCPIP titration solution at pH 8', display_color="#00A6503E")

    for well in buffer_pH6:
        well.load_liquid(liquid=buffer_pH6_liquid, volume=15000)
    for well in buffer_pH7:
        well.load_liquid(liquid=buffer_pH7_liquid, volume=15000)
    for well in buffer_pH8:
        well.load_liquid(liquid=buffer_pH8_liquid, volume=15000)
    for well in dcpip_pH6:
        well.load_liquid(liquid=dcpip_pH6_liquid, volume=300)
    for well in dcpip_pH7:
        well.load_liquid(liquid=dcpip_pH7_liquid, volume=300)
    for well in dcpip_pH8:
        well.load_liquid(liquid=dcpip_pH8_liquid, volume=300)

def add_20uL_buffer(protocol):
    """Add the initial buffer volume to every non-DCPIP assay column.

    DCPIP-only wells in columns 1, 9, and 17 receive no buffer. The final buffer-only
    wells in columns 8, 16, and 24 start with 40 uL to keep the final 60 uL total
    consistent with the titration design.
    """
    # pH 6 block: columns 1-8 (column 1 has no buffer; column 8 gets 40 uL)
    p300m.pick_up_tip()
    for idx, column in enumerate(plate.columns()[0:8]):
        if idx == 0:
            continue
        volume = 40 if idx == 7 else 10
        p300m.distribute(volume, buffer_pH6[0], column, new_tip='never')
    p300m.drop_tip()

    # pH 7 block: columns 9-16 (column 9 has no buffer; column 16 gets 40 uL)
    p300m.pick_up_tip()
    for idx, column in enumerate(plate.columns()[8:16]):
        if idx == 0:
            continue
        volume = 40 if idx == 7 else 10
        p300m.distribute(volume, buffer_pH7[0], column, new_tip='never')
    p300m.drop_tip()

    # pH 8 block: columns 17-24 (column 17 has no buffer; column 24 gets 40 uL)
    p300m.pick_up_tip()
    for idx, column in enumerate(plate.columns()[16:24]):
        if idx == 0:
            continue
        volume = 40 if idx == 7 else 10
        p300m.distribute(volume, buffer_pH8[0], column, new_tip='never')
    p300m.drop_tip()


def add_and_titrate_dcpip(protocol):
    columns_pH6 = plate.columns()[1:7]
    columns_pH7 = plate.columns()[9:15]
    columns_pH8 = plate.columns()[17:23]
    dcpip_conditions = (
        (dcpip_pH6, plate.columns()[0], columns_pH6),
        (dcpip_pH7, plate.columns()[8], columns_pH7),
        (dcpip_pH8, plate.columns()[16], columns_pH8),
    )

    # Add DCPIP once and then serially transfer within that pH condition using the same tip.
    for dcpip_source, dcpip_only_column, titration_columns in dcpip_conditions:
        p300m.pick_up_tip()
        p300m.distribute(90, dcpip_source[0], dcpip_only_column, new_tip='never')
        serial_sources = [dcpip_only_column, *titration_columns[:-1]]
        p300m.transfer(30, serial_sources, titration_columns,
                       new_tip='never', mix_after=(3, 30))
        p300m.drop_tip()

    # Add the final 50 uL to every assay well except the DCPIP-only wells,
    # then remove the extra 10 uL from the final titration column in each pH block.
    buffer_groups = (
        (buffer_pH6[0], plate.columns()[1:8], plate.columns()[6]),
        (buffer_pH7[0], plate.columns()[9:16], plate.columns()[14]),
        (buffer_pH8[0], plate.columns()[17:24], plate.columns()[22]),
    )
    for buffer_source, columns, excess_column in buffer_groups:
        p300m.pick_up_tip()
        for column in columns:
            if column == plate.columns()[0] or column == plate.columns()[8] or column == plate.columns()[16]:
                continue
            p300m.distribute(50, buffer_source, column, new_tip='never', mix_after=(3, 50))
        p300m.aspirate(10, excess_column[0])
        p300m.dispense(10, liquid_waste)
        p300m.drop_tip()

