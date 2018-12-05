from opentrons import robot,labware, instruments
import math

###This has been tested and worked well with 10 ul tips, when pipetting with 10 ul tip onlyself.
##Do not use two pipettes at the same time at the moment to avoid tips crashing the plateself.


example_csv = """
A1, 20
A3, 10
B2, 15

"""

ep_rack = labware.load('epMotion2', '1')

tiprack_10 = [labware.load('tiprack-10ul', slot)
              for slot in ['6', '7']]
tiprack_50 = labware.load('tiprack-200ul', '9')


pool=ep_rack.cols('1')

p10 = instruments.P10_Single(
    mount='left',
    tip_racks=tiprack_10)

p50 = instruments.P50_Single(
    mount='right',
    tip_racks=tiprack_50)

######### Pooling
#1.Get the number and location from CSV files
#2.Transfer a particular volume from particular well to the epMotion reservoir 10 mL
#3.Change tips
#The current dead volume for PCR-flat is 3 ul.
#For the small volume, do one in two dilutions in those. Then, double the amount of samples taken. => do manually.

#Do one plate at a time, CSV one plate at a time, change tips. Then move to another plate.

def run_custom_protocol(
        volumes_csv: 'FileInput'=example_csv,
        tip_reuse: 'StringSelection...'='new tip each time'
        ):
    data = [
        [well, float(vol)]
        for well, vol in
        [row.split(',') for row in volumes_csv.strip().split('\n') if row]
    ]

#    source_plate = sample_plate_1
    dest_plate = pool

#    tip_strategy = 'always' if tip_reuse == 'new tip each time' else 'once'
    for well_idx, (source_well, vol) in enumerate(data):
        if vol>10:
            p50.transfer(
            vol,
            source_plate.wells(source_well),
            pool,
            new_tip='always')
        else:
            p10.transfer(
            vol,
            source_plate.wells(source_well),
            pool,
            new_tip='always')

##SS Plate 1
source_plate = labware.load('96-flat','3')
run_custom_protocol(**{'volumes_csv': 'A12,0.794\r\nB1,0.938\r\nD4,0.966\r\nD2,1.007\r\nC11,1.013\r\nH1,1.034\r\nF8,1.057\r\nH5,1.059\r\nC5,1.066\r\nD3,1.078\r\nD5,1.088\r\nF4,1.095\r\nB12,1.098\r\nF5,1.111\r\nB8,1.112\r\nH11,1.132\r\nD11,1.160\r\nF6,1.161\r\nB5,1.174\r\nD9,1.199\r\nF11,1.214\r\nD8,1.217\r\nA6,1.218\r\nD1,1.228\r\nH6,1.229\r\nA5,1.236\r\nC8,1.239\r\nB4,1.249\r\nD6,1.260\r\nF7,1.265\r\nF3,1.271\r\nC6,1.288\r\nH2,1.289\r\nB6,1.299\r\nH12,1.303\r\nB11,1.304\r\nB10,1.311\r\nB2,1.316\r\nB7,1.323\r\nA4,1.328\r\nC2,1.330\r\nC4,1.336\r\nA7,1.365\r\nD10,1.391\r\nA1,1.399\r\nH9,1.406\r\nC7,1.409\r\nF2,1.417\r\nC10,1.464\r\nA3,1.482\r\nH8,1.482\r\nC12,1.486\r\nD7,1.506\r\nA9,1.506\r\nF9,1.507\r\nH10,1.507\r\nA10,1.534\r\nF12,1.536\r\nB9,1.567\r\nD12,1.573\r\nA8,1.593\r\nA11,1.604\r\nG3,1.625\r\nC3,1.697\r\nG8,1.723\r\nH7,1.737\r\nG4,1.817\r\nG5,1.822\r\nB3,1.844\r\nE7,1.872\r\nA2,1.916\r\nG7,1.971\r\nF1,2.016\r\nF10,2.043\r\nG9,2.286\r\nC1,2.412\r\nG10,2.686\r\nE5,2.901\r\nE11,3.388\r\nG2,4.046\r\nE3,4.141\r\nE1,4.227\r\nE2,4.346\r\nE6,5.642\r\nE10,5.802\r\nE4,6.578\r\nE12,6.758\r\nG12,7.069\r\nE8,7.593\r\nG6,7.935\r\nG1,8.662\r\nE9,8.913\r\nG11,9.179', 'tip_reuse': 'reuse tip'})


##SS Plate 2
source_plate = labware.load('96-flat','4')
run_custom_protocol(**{'volumes_csv': 'A12,0.967\r\nH12,1.041\r\nD1,1.043\r\nG11,1.046\r\nC5,1.059\r\nH3,1.060\r\nD4,1.064\r\nG12,1.067\r\nG1,1.073\r\nC7,1.075\r\nF6,1.076\r\nH1,1.091\r\nH10,1.091\r\nE9,1.103\r\nG3,1.114\r\nE10,1.114\r\nF7,1.121\r\nD9,1.133\r\nG6,1.150\r\nF4,1.151\r\nD5,1.151\r\nB9,1.157\r\nB11,1.160\r\nE3,1.164\r\nD12,1.170\r\nB1,1.180\r\nE1,1.180\r\nE6,1.191\r\nD11,1.196\r\nE4,1.198\r\nH4,1.202\r\nF11,1.204\r\nG2,1.213\r\nB12,1.214\r\nG4,1.230\r\nE12,1.236\r\nE5,1.244\r\nB7,1.246\r\nB8,1.249\r\nC4,1.253\r\nF12,1.263\r\nA8,1.265\r\nD6,1.271\r\nH11,1.271\r\nC6,1.275\r\nC11,1.276\r\nA3,1.280\r\nA6,1.281\r\nE8,1.283\r\nC9,1.302\r\nB4,1.317\r\nD10,1.320\r\nF1,1.323\r\nA1,1.330\r\nC3,1.336\r\nH9,1.336\r\nA5,1.337\r\nB5,1.337\r\nC10,1.338\r\nG5,1.341\r\nD7,1.349\r\nB6,1.350\r\nF3,1.352\r\nC8,1.352\r\nA2,1.358\r\nA7,1.387\r\nH2,1.388\r\nG10,1.390\r\nB2,1.411\r\nH5,1.434\r\nG9,1.440\r\nB10,1.447\r\nA10,1.450\r\nF9,1.459\r\nF5,1.461\r\nE11,1.461\r\nF10,1.470\r\nC12,1.470\r\nG8,1.507\r\nD8,1.515\r\nA9,1.519\r\nG7,1.538\r\nH8,1.608\r\nE7,1.621\r\nC2,1.690\r\nH6,1.735\r\nF8,1.775\r\nA11,1.806\r\nE2,2.139\r\nC1,2.189\r\nF2,2.398\r\nB3,2.563\r\nD2,6.440\r\nD3,5.834', 'tip_reuse': 'reuse tip'})

##SS Plate 3
#run_custom_protocol(**{'sample_plate_3','volumes_csv': 'A2,0.510,1.0\r\nD4,0.682,1.4\r\nA1,0.885,1.8\r\nA6,0.941,1.9\r\nH1,0.947,1.9\r\nF7,1.047,2.1\r\nE7,1.101,2.2\r\nG2,1.126,2.3\r\nB7,1.146,2.3\r\nC5,1.170,2.3\r\nA4,1.234,2.5\r\nF5,1.246,2.5\r\nG7,1.289,2.6\r\nA5,1.295,2.6\r\nF6,1.317,2.6\r\nC2,1.340,2.7\r\nE4,1.365,2.7\r\nB1,1.393,2.8\r\nD7,1.412,2.8\r\nF1,1.422,2.8\r\nC4,1.424,2.8\r\nB4,1.435,2.9\r\nE1,1.477,3.0\r\nD5,1.573,3.1\r\nD1,1.597,3.2\r\nE5,1.647,3.3\r\nE2,1.664,3.3\r\nE6,1.723,3.4\r\nF4,1.777,3.6\r\nB5,1.801,3.6\r\nC7,1.817,3.6\r\nG6,1.841,3.7\r\nA3,1.861,3.7\r\nC3,1.866,3.7\r\nF2,1.884,3.8\r\nG3,1.925,3.8\r\nC1,1.994,4.0\r\nF3,1.994,4.0\r\nD3,2.033,4.1\r\nH2,2.121,4.2\r\nD2,2.189,4.4\r\nB2,2.253,4.5\r\nC6,2.334,4.7\r\nA7,2.343,4.7\r\nH5,2.795,5.6\r\nH3,2.854,5.7\r\nG4,2.922,5.8\r\nB6,2.943,5.9\r\nD6,2.964,5.9\r\nH6,3.333,6.7\r\nH4,3.455,6.9\r\nG5,3.565,7.1\r\nH7,3.639,7.3\r\nE3,3.739,7.5\r\nG1,4.980,10.0\r\nB3,5.775,11.5', 'pipette_axis': 'B (left side)', 'pipette_model': 'p1000', 'source_plate_type': '96-flat', 'destination_plate_type': '96-flat', 'tip_reuse': 'reuse tip'})

robot.home()
