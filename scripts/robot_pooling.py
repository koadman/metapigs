#!/usr/bin/env python3

import numpy as np
import sys
import re
#from opentrons import containers, instruments

# alter the ratio of the following samples:
custom={
    "7D9":10,
    "8G12":10,
    "5A2":10
}

# configure the unused well blacklist
blacklist={}
for col in range(1,12):
    blacklist['2D'+str(col)]=1

for col in range(9,12):
    blacklist['2H'+str(col)]=1

for plate in [9,10]:
    for col in [11,12]:
        for row in ['A','B','C','D','E','F','G','H']:
            blacklist[str(plate)+row+str(col)]=1


# some parameters for pipetting
max_pool_vol = 250
aliquot_vol = 10
samples_per_pool = max_pool_vol / aliquot_vol

rc_floor = 50 # treat all samples as though they have this many reads
manual_min = 2.5 # minimum amount we would manually pipet

# read the count data
count_table = np.genfromtxt(sys.argv[1],skip_header=1,dtype=str)

# convert names to plate codes
well_counts = {}
for s_row in range(count_table.shape[0]):
    sample = count_table[s_row,0]
    rrr = re.search('plate_(\d+)_(\w)(\d+)', sample)
    plate = rrr.group(1)
    row = rrr.group(2)
    col = rrr.group(3)
    code = str(plate)+row+str(col)
    well_counts[code]=max(rc_floor,int(count_table[s_row,1]))
    if code in custom:
        well_counts[code] /= custom[code]

# add zero counts for any missing wells
for plate in range(1,10):
    for row in ['A','B','C','D','E','F','G','H']:
        for col in range(1,12):
            code = str(plate)+row+str(col)
            if not code in well_counts:
                well_counts[code]=rc_floor

cur_pool = 0
sample_count = 0
pool_sample_count = {cur_pool:0}
pool_read_count = {cur_pool:0}
pool_minmax = {cur_pool:[rc_floor,rc_floor]}
max_relrange = 1.1
pool_samples = {0:[]}
for count, well in sorted((value, key) for (key,value) in well_counts.items()):
    if ( pool_sample_count[cur_pool] >= samples_per_pool or
        (cur_pool > 5 and pool_minmax[cur_pool][0]*max_relrange < count) or
        (pool_minmax[cur_pool][0] < 140 and count > 140) or
        (pool_minmax[cur_pool][0] < 250 and count > 250)):
        # start a new pool
        cur_pool += 1
        pool_sample_count[cur_pool] = 0
        pool_read_count[cur_pool] = 0
        pool_minmax[cur_pool] = [count,0]
        pool_samples = {cur_pool:[]}
    if well in blacklist:
        if count > 0:
            print("# WARNING: blacklisted well ID "+well+" has read count "+str(count))
    else:
        # add a sample to the pool
        pool_sample_count[cur_pool] += 1
        pool_read_count[cur_pool] += count
        pool_minmax[cur_pool][1] = count
        pool_samples[cur_pool].append(well)
        print("well ID "+well+" read count "+str(count)+" into pool "+str(cur_pool))

# pool 0 is least concentrated - use all of this one and less of others
pool_0_ul = pool_sample_count[0]*aliquot_vol
pool_0_avgconc = pool_read_count[0] / pool_sample_count[0]
for pool in pool_sample_count:
    pool_avgconc = pool_read_count[pool] / pool_sample_count[pool]
#    print(str(pool_0_avgconc)+"\t"+str(pool_minmax[pool][1]))
    pool_relconc = pool_0_avgconc / pool_minmax[pool][1]
    pool_ul = pool_0_ul * pool_relconc
    print("Aliquot "+str(pool_ul)+"ul from pool "+str(pool)+" poolvol " + str(pool_sample_count[pool]*10) + " relrange "+
        str(pool_minmax[pool][0]/pool_minmax[pool][1]))


# robot instructions start here


p10rack = containers.load('tiprack-10ul', 'B2', 'p10rack')

# look here for dimensions https://docs.opentrons.com/ot1/containers.html
lib_plates = {
	1:containers.load('96-PCR-tall', 'A1', 'lib_plate_1'),
	2:containers.load('96-PCR-tall', 'A2', 'lib_plate_2'),
	3:containers.load('96-PCR-tall', 'A3', 'lib_plate_3'),
	4:containers.load('96-PCR-tall', 'B1', 'lib_plate_4'),
	5:containers.load('96-PCR-tall', 'B3', 'lib_plate_5'),
	6:containers.load('96-PCR-tall', 'C1', 'lib_plate_6'),
	7:containers.load('96-PCR-tall', 'C3', 'lib_plate_7'),
	8:containers.load('96-PCR-tall', 'D1', 'lib_plate_8'),
	9:containers.load('96-PCR-tall', 'D2', 'lib_plate_9'),
	10:containers.load('96-PCR-tall', 'D3', 'lib_plate_10')
}

pool_plate = containers.load('96-PCR-flat', 'C2', 'pool_plate')
trash = containers.load('trash-box', 'E1', 'trash')

cur_row=1
row_letters={1:'A',2:'B',3:'C',4:'D',5:'E',6:'F',7:'G',8:'H'}
cur_col=1
tip_col=1
tip_row=1
for pool in pool_samples:
    pool_from = []
    for sample in pool_samples[pool]:
        rrr = re.search('(\d+)(\w\d+)',sample)
        plate = rrr.group(1)
        well = rrr.group(2)
        pool_from.append( lib_plates[ plate ].well(well).bottom() )
    pool_dest_well = row_letters[cur_row]+str(cur_col)
    pool_dest = pool_plate.well(pool_dest_well).bottom()
    cur_col += 1
    if cur_col > 12:
        cur_col = 1
        cur_row += 1

    #pick up new tip first
    tip_well = row_letters[tip_row]+str(tip_col)
    tip_col += 1
    if tip_col > 12:
        tip_col = 1
        tip_row += 1

    pipette.pick_up_tip(p10rack.wells(tip_well))
    p10.transfer(
        10,
        pool_from,
        pool_dest,
        disposal_vol=0,
        mix_before=(3), # mix 3 times
        mix_after=(3),  # mix 3 times
        touch_tip=True,
        blow_out=True,
        new_tip='never')

    pipette.drop_tip(trash)
    robot.home()
