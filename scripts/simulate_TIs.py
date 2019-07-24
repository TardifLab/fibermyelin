#!/usr/bin/python
#/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tues Jun 26 2019

@author: jcampbel
"""

import numpy as np
import nibabel as nib
import argparse
import math
import subprocess
import matplotlib.pyplot as plt
import os
import argparse



#run from /data/mril/mril4/jcampbel/IR-diff_test/TI_sim

#make a mask with singles and crossings, ~100 voxels
#in fsleyes:
#Overlay-copy-mask w same dims
#swtich to copy - edit mode-
#create mask...

#fslmaths index_firstframe_strides.nii.gz -thr 1.5 -mas ../irdiff_4v8_20190319_121513/mb3/masks/mask3.nii.gz mask3_crossings
#fslmaths mask3_crossings.nii.gz -bin mask3
#fslmaths index_firstframe_strides.nii.gz  -thr 0.5 -uthr 1.5 -mas ../irdiff_4v8_20190319_121513/mb3/masks/mask2.nii.gz mask2_singles.nii
#fslmaths mask3_singles.nii.gz -bin mask3_singles_bin.nii.gz

global crossings_different_T1s
crossings_different_T1s=True

parser = argparse.ArgumentParser(description='')

parser.add_argument('-o', dest='myresult_output_filename', help='output filename (.npy file). Optional', required=False, nargs=1)
parser.add_argument('-load', dest='myresult_input_filename', help='input filename of previously computed result, to plot it. Optional', required=False, nargs=1)

myargs = parser.parse_args()



if (myargs.myresult_input_filename):
    myresult = np.load(myargs.myresult_input_filename)
    print myresult

else:

    num_shuffles = 3

    num_SMS = 4
    min_TI = 25
    max_TI = 3000
    num_TIs = num_SMS*num_shuffles
    num_slices = 72

    gap_step=300

    num_gaps = int((int(math.floor(max_TI / 12))-110)/gap_step) + 1

    myresult=np.zeros([num_gaps, num_shuffles, 2, 4])

    for gap_index in range(num_gaps):
        gap = 110 + gap_index * gap_step
        TIs = np.zeros([num_shuffles, num_TIs / num_shuffles])
        for shuffle in range (0, num_shuffles):
            f= open("TIsim.txt","w+")
            for slice  in range (0, num_slices):
                for i in range(int(num_TIs / num_shuffles)): # has to be int!
                    TI = min_TI + (shuffle + num_shuffles * i) * gap
                    f.write("%f " % TI)
                f.write("\n")
            f.close()


            #run the sim:

            subprocess.call(['/home/bic/jcampbel/source/fibermyelin_jc/fibermyelin/scripts/fibermyelin_pipeline3.py',
                                '-t1', 't1_test',
                                '-Dpar', 'Dpar_test.nii',
                                '-mask', '/data/mril/mril4/jcampbel/IR-diff_test/irdiff_4v8_20190319_121513/mb3/masks/mask3.nii.gz',
                                '-vic', '/data/mril/mril4/jcampbel/IR-diff_test/irdiff_4v8_20190319_121513/mb3/hardi-analysis/AMICO/NODDI/FIT_ICVF-strides.nii.gz',
                                '-viso', '/data/mril/mril4/jcampbel/IR-diff_test/irdiff_4v8_20190319_121513/mb3/FIT_ISOVF_strides_xflip.nii',
                                '-afd', '/data/mril/mril4/jcampbel/IR-diff_test/irdiff_4v8_20190319_121513/mb3/hardi-analysis/afd_voxel_strides.nii',
                                '-afdthresh', '0.2',
                                '-dirs', '/data/mril/mril4/jcampbel/IR-diff_test/irdiff_4v8_20190319_121513/mb3/hardi-analysis/directions_voxel_strides.nii',
                                '-IRdiff', '/data/mril/mril4/jcampbel/IR-diff_test/irdiff_4v8_20190319_121513/mb3/ir-analysis/ir-diff-strides.nii',
                                '-TIs', 'TIsim.txt',
                                '-bvals', '/data/mril/mril4/jcampbel/IR-diff_test/irdiff_4v8_20190319_121513/mb3/bval_rmdummy.bval',
                                '-bvecs', '/data/mril/mril4/jcampbel/IR-diff_test/irdiff_4v8_20190319_121513/mb3/bvec_rmdummy.bvec',
                                '-TR', '3000',
                                '-TE', '56',
                                '-fixel', 'fixel_test'])


            #break into singles and crossings:
            subprocess.call(['mrcalc', 't1_test.nii', '/data/mril/mril4/jcampbel/IR-diff_test/TI_sim/mask3_singles_bin.nii.gz', '-mult', 't1_test_singles.nii', '-force'])
            subprocess.call(['mrcalc', 't1_test.nii', '/data/mril/mril4/jcampbel/IR-diff_test/TI_sim/mask3_crossings_bin.nii.gz', '-mult', 't1_test_crossings.nii', '-force'])

            #separate the first and second fiber:
            subprocess.call(['fslroi', 't1_test_crossings.nii', 't1_test_crossings_firstfiber.nii','0', '1'])
            subprocess.call(['fslroi', 't1_test_crossings.nii', 't1_test_crossings_secondfiber.nii','1', '1'])

            #get mean and std from fslstats, for singles and crossings:


            myresult[gap_index, shuffle, 0, 0] = subprocess.check_output(['fslstats', 't1_test_singles.nii', '-M'])
            myresult[gap_index, shuffle, 0, 1] = subprocess.check_output(['fslstats', 't1_test_crossings.nii', '-M'])
            myresult[gap_index, shuffle, 0, 2] = subprocess.check_output(['fslstats', 't1_test_crossings_firstfiber.nii', '-M'])
            myresult[gap_index, shuffle, 0, 3] = subprocess.check_output(['fslstats', 't1_test_crossings_secondfiber.nii', '-M'])
            myresult[gap_index, shuffle, 1, 0] = subprocess.check_output(['fslstats', 't1_test_singles.nii', '-S'])
            myresult[gap_index, shuffle, 1, 1] = subprocess.check_output(['fslstats', 't1_test_crossings.nii', '-S'])
            myresult[gap_index, shuffle, 1, 2] = subprocess.check_output(['fslstats', 't1_test_crossings_firstfiber.nii', '-S'])
            myresult[gap_index, shuffle, 1, 3] = subprocess.check_output(['fslstats', 't1_test_crossings_secondfiber.nii', '-S'])


    np.save(myargs.myresult_output_filename[0], np.array(myresult))

print myresult
#plot:

#bias in singles
plt.figure(1, figsize=(9, 6))
plt.subplot(241)

#DO: 110 + gap_step * range(num_gaps) fails, currently not in ms, just a gap index

for shuffle in range (0, num_shuffles):
    plt.plot(range(num_gaps), myresult[:, shuffle, 0, 0]-750,label=['shuffle %d' % shuffle])
plt.legend(loc='best')
plt.title('Bias; single fibers', fontsize=18)
plt.xlabel('gap (ms)')
plt.ylabel('bias (ms)')

if not crossings_different_T1s:
    #bias in crossings
    plt.subplot(242)
    for shuffle in range (0, num_shuffles):
        plt.plot(range(num_gaps), myresult[:, shuffle, 0, 1]-750,label=['shuffle %d' % shuffle])

    plt.legend(loc='best')
    plt.title('Bias; crossing fiber', fontsize=18)
    plt.xlabel('gap (ms)')
    plt.ylabel('bias (ms)')

else:
    #bias in first fiber
    plt.subplot(243)
    for shuffle in range (0, num_shuffles):
        plt.plot(range(num_gaps), myresult[:, shuffle, 0, 2]-700,label=['shuffle %d' % shuffle])

    plt.legend(loc='best')
    plt.title('Bias; first fiber', fontsize=18)
    plt.xlabel('gap (ms)')
    plt.ylabel('bias (ms)')


    #bias in second fiber
    plt.subplot(244)
    for shuffle in range (0, num_shuffles):
        plt.plot(range(num_gaps), myresult[:, shuffle, 0, 3]-800,label=['shuffle %d' % shuffle])

    plt.legend(loc='best')
    plt.title('Bias; second fiber', fontsize=18)
    plt.xlabel('gap (ms)')
    plt.ylabel('bias (ms)')

#std in singles
plt.subplot(245)
for shuffle in range (0, num_shuffles):
    print
    plt.plot(range(num_gaps),myresult[:, shuffle, 1, 0],label=['shuffle %d' % shuffle])
plt.legend(loc='best')
plt.title('Standard deviation; single fibers', fontsize=18)
plt.xlabel('gap (ms)')
plt.ylabel('standard deviation (ms)')

if not crossings_different_T1s:
    #std in crossings
    plt.subplot(246)
    for shuffle in range (0, num_shuffles):
        plt.plot(range(num_gaps),myresult[:, shuffle, 1, 1],label=['shuffle %d' % shuffle])
    plt.legend(loc='best')
    plt.title('Standard deviation; crossing fibers', fontsize=18)
    plt.xlabel('gap (ms)')
    plt.ylabel('standard deviation (ms)')

else:
    #std in first fiber
    plt.subplot(247)
    for shuffle in range (0, num_shuffles):
        plt.plot(range(num_gaps),myresult[:, shuffle, 1, 2],label=['shuffle %d' % shuffle])
    plt.legend(loc='best')
    plt.title('Standard deviation; first fiber', fontsize=18)
    plt.xlabel('gap (ms)')
    plt.ylabel('standard deviation (ms)')

    #std in second fiber
    plt.subplot(248)
    for shuffle in range (0, num_shuffles):
        plt.plot(range(num_gaps),myresult[:, shuffle, 1, 3],label=['shuffle %d' % shuffle])
    plt.legend(loc='best')
    plt.title('Standard deviation; second fiber', fontsize=18)
    plt.xlabel('gap (ms)')
    plt.ylabel('standard deviation (ms)')

plt.show()

plt.savefig("sim_result_fig.png")