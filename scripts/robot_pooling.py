#!/usr/bin/env python3

import sys
import re
from opentrons import robot, containers, instruments

robot.get_warnings()

# alter the ratio of the following samples:
custom={
    "7D9":10,
    "8G12":10,
    "5A2":10
}

# configure the unused well blacklist
blacklist={}
for col in range(1,13):
    blacklist['2D'+str(col)]=1

for col in range(9,13):
    blacklist['2H'+str(col)]=1

for plate in [9,10]:
    for col in [11,12]:
        for row in ['A','B','C','D','E','F','G','H']:
            blacklist[str(plate)+row+str(col)]=1


# samples which lack data for normalisation go to a separate pool

awol_pool={'1A1','2A1','2B1','2C1','6A7','4F10','8E11'}


# some parameters for pipetting
max_pool_vol = 250
aliquot_vol = 10
samples_per_pool = max_pool_vol / aliquot_vol

rc_floor = 50 # treat all samples as though they have this many reads
manual_min = 2.5 # minimum amount we would manually pipet


#replaced count_file = open(sys.argv[1]) with count_file = open('/Users/12705859/metapigs/source_data/plate_counts.tsv')

# read the count data
count_file = open('/Users/12705859/metapigs/source_data/plate_counts.tsv')
count_table = {}
for line in count_file:
    ll = line.rstrip().split()
    count_table[ll[0]]=int(ll[1])

# convert names to plate codes
well_counts = {}
for sample in count_table:
    rrr = re.search('plate_(\d+)_(\w)(\d+)', sample)
    plate = rrr.group(1)
    row = rrr.group(2)
    col = rrr.group(3)
    code = str(plate)+row+str(col)
    if code not in blacklist:
        well_counts[code]=max(rc_floor,count_table[sample])
    else:
        well_counts[code]=count_table[sample]

    if code in custom:
        well_counts[code] /= custom[code]

# add zero counts for any missing wells
for plate in range(1,10):
    for row in ['A','B','C','D','E','F','G','H']:
        for col in range(1,12):
            code = str(plate)+row+str(col)
            if not code in well_counts:
                if not code in blacklist:
                    print("sample "+code+" was missing a read count")
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
        pool_samples[cur_pool]=[]
    if (well not in blacklist) and (well not in awol_pool):
        # add a sample to the pool
        pool_sample_count[cur_pool] += 1
        pool_read_count[cur_pool] += count
        pool_minmax[cur_pool][1] = count
        pool_samples[cur_pool].append(well)
        print("well ID "+well+" read count "+str(count)+" into pool "+str(cur_pool))

# make the awol pool
cur_pool += 1
pool_samples[cur_pool] = list(awol_pool)
pool_sample_count[cur_pool] = len(list(awol_pool))
pool_read_count[cur_pool] = 3000 * pool_sample_count[cur_pool]
pool_minmax[cur_pool] = [3000,3000]

# pool 0 is least concentrated - use all of this one and less of others
pool_0_ul = pool_sample_count[0]*aliquot_vol
pool_0_avgreads = pool_read_count[0] / pool_sample_count[0]
water_pools = {}
for pool in pool_sample_count:
    pool_avgreads = pool_read_count[pool] / pool_sample_count[pool]
    pool_relconc = pool_0_avgreads / pool_avgreads
    pool_ul = pool_sample_count[pool] * aliquot_vol * pool_relconc
    # dilute with water if we would do a tiny aliquot and the well has room
    if(pool_ul < 3 and pool_sample_count[pool] < samples_per_pool):
        pool_ul = pool_sample_count[pool] * aliquot_vol * pool_0_avgreads / (pool_read_count[pool] / samples_per_pool)
        water_pools[pool]='yes'

    print("Aliquot "+str(pool_ul)+"ul from pool "+str(pool)+" samples " + str(pool_sample_count[pool]) + " avgreads " + str(pool_avgreads) + " relrange "+
        str(pool_minmax[pool][0]/pool_minmax[pool][1]))

#
# robot instructions start here
#
#p10 = containers.load('tiprack-10ul', 'B2', 'p10rack') changed into:
p10rack = containers.load('tiprack-10ul', 'D2', 'p10rack')



trash_container = containers.load('trash-box', 'E1', 'trash')

#added this description of pipette
p10 = instruments.Pipette(
    axis='a',
    name='p10',
    max_volume=10,
    min_volume=0.5,
    channels=1,
    trash_container=trash_container,
    tip_racks=[p10rack])

# look here for dimensions https://docs.opentrons.com/ot1/containers.html
lib_plates = {
	1:containers.load('96-PCR-tall', 'B1', 'lib_plate_1'),
#	2:containers.load('96-PCR-tall', 'B2', 'lib_plate_2'),
#	3:containers.load('96-PCR-tall', 'B3', 'lib_plate_3'),
#	4:containers.load('96-PCR-tall', 'C3', 'lib_plate_4'),
#	5:containers.load('96-PCR-tall', 'D3', 'lib_plate_5'),
#	6:containers.load('96-PCR-tall', 'B1', 'lib_plate_6'),
#	7:containers.load('96-PCR-tall', 'B2', 'lib_plate_7'),
#	8:containers.load('96-PCR-tall', 'B3', 'lib_plate_8'),
#	9:containers.load('96-PCR-tall', 'C3', 'lib_plate_9'),
#	10:containers.load('96-PCR-tall', 'D3', 'lib_plate_10')
}

pool_plate = containers.load('96-PCR-flat', 'C2', 'pool_plate')

water = containers.load('trough-12row', 'C1')
water_well = water.wells('A1')

# preload some wells with water before the first batch of plates is processed
# so that high concentration pools don't require tiny aliquots
cur_row=1
cur_col=1
row_letters={1:'A',2:'B',3:'C',4:'D',5:'E',6:'F',7:'G',8:'H'}
if 1 in lib_plates:
    for pool in pool_samples:
        if pool in water_pools:
            pool_dest_well = row_letters[cur_row]+str(cur_col)
            pool_dest = pool_plate.well(pool_dest_well)
            # pipette 10uL of water for each sample less than the samples per pool
            print("pipetting "+str(int(samples_per_pool - pool_sample_count[pool])*10)+"ul water for pool "+str(pool)+" well "+pool_dest_well)
            for i in range(int(samples_per_pool - pool_sample_count[pool])):
                p10.transfer(
                    10,
                    water_well,
                    pool_dest,
                    disposal_vol=0,
                    mix_before=(0),
                    mix_after=(0),
                    touch_tip=True,
                    blow_out=True,
                    new_tip='never')
        # advance to next pool
        cur_col += 1
        if cur_col > 12:
            cur_col = 1
            cur_row += 1

# start making the pools
cur_row=1
cur_col=1
for pool in pool_samples:
    pool_from = []
    for sample in pool_samples[pool]:
        rrr = re.search('(\d+)(\w\d+)',sample)
        plate = int(rrr.group(1))
        well = rrr.group(2)
        if plate in lib_plates:
            pool_from.append( lib_plates[ plate ].well(well) )
    pool_dest_well = row_letters[cur_row]+str(cur_col)
    pool_dest = pool_plate.well(pool_dest_well)
    cur_col += 1
    if cur_col > 12:
        cur_col = 1
        cur_row += 1
    # ignore this pool if empty
    if len(pool_from) == 0:
        continue
        
    p10.transfer(
        10,
        pool_from,
        pool_dest,
        disposal_vol=0,
        mix_before=(2, 10), # mix 3 times
        mix_after=(3),  # mix 3 times
        touch_tip=True,
        blow_out=True,
        new_tip='never')

for c in robot.commands():
    print(c)
        