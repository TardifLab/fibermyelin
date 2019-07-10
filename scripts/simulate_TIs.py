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

    gap_step=10

    num_gaps = int((int(math.floor(max_TI / 12))-110)/gap_step) + 1

    myresult=np.zeros([num_gaps, num_shuffles, 2, 2])

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

            subprocess.call(['/home/bic/jcampbel/source/fibermyelin_jc/fibermyelin/scripts/fibermyelin_pipeline.py',
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
            subprocess.call(['fslmaths', 't1_test.nii', '-mul', '/data/mril/mril4/jcampbel/IR-diff_test/TI_sim/mask3_singles_bin.nii.gz', 't1_test_singles.nii'])
            subprocess.call(['fslmaths', 't1_test.nii','-mul', '/data/mril/mril4/jcampbel/IR-diff_test/TI_sim/mask3_crossings_bin.nii.gz', 't1_test_crossings.nii'])



            #get mean and std from fslstats,for singles and crossings:


            myresult[gap_index, shuffle, 0, 0] = subprocess.check_output(['fslstats', 't1_test_singles.nii.gz', '-M'])
            myresult[gap_index, shuffle, 0, 1] = subprocess.check_output(['fslstats', 't1_test_crossings.nii.gz', '-M'])
            myresult[gap_index, shuffle, 1, 0] = subprocess.check_output(['fslstats', 't1_test_singles.nii.gz', '-S'])
            myresult[gap_index, shuffle, 1, 1] = subprocess.check_output(['fslstats', 't1_test_crossings.nii.gz', '-S'])


    np.save("myargs.myresult_output_filename", np.array(myresult))

print myresult
#plot:

#bias in singles
plt.figure(1, figsize=(9, 6))
plt.subplot(221)

#DO: 110 * gap_step * range(num_gaps) fails, currently not in ms, just a gap index

for shuffle in range (0, num_shuffles):
    plt.plot(range(num_gaps), myresult[:, shuffle, 0, 0]-750,label=['shuffle %d' % shuffle])
plt.legend(loc='best')
plt.title('Bias; single fibres', fontsize=18)
plt.xlabel('gap (ms)')
plt.ylabel('bias (s)')

#bias in crossings
plt.subplot(222)
for shuffle in range (0, num_shuffles):
    plt.plot(range(num_gaps), myresult[:, shuffle, 0, 1]-750,label=['shuffle %d' % shuffle])
    #ax.plot([110 * np.ones(num_gaps) + range(num_gaps) * 25], myresult[:, shuffle, 0, 1]-750,label=['shuffle %d' % shuffle])
plt.legend(loc='best')
plt.title('Bias; crossing fibres', fontsize=18)
plt.xlabel('gap (ms)')
plt.ylabel('bias (s)')


#bias in singles
plt.subplot(223)
for shuffle in range (0, num_shuffles):
    print
    plt.plot(range(num_gaps),myresult[:, shuffle, 1, 0],label=['shuffle %d' % shuffle])
plt.legend(loc='best')
plt.title('Standard deviation; single fibres', fontsize=18)
plt.xlabel('gap (ms)')
plt.ylabel('standard deviation (s)')

#bias in crossings
plt.subplot(224)
for shuffle in range (0, num_shuffles):
    plt.plot(range(num_gaps),myresult[:, shuffle, 1, 1],label=['shuffle %d' % shuffle])
plt.legend(loc='best')
plt.title('Standard deviation; crossing fibres', fontsize=18)
plt.xlabel('gap (ms)')
plt.ylabel('standard deviation (s)')



plt.show()

