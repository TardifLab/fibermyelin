#!/usr/bin/python
#/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tues Jun 26 2019
Modified on Feb 25th 2020 for MTRs
 (@ileppe)
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



#run from /data_/tardiflab/ilana/mt-diff/sim
# set this is FiberMTSolver.py simulate=True
#make a mask with singles and crossings, ~100 voxels
#in fsleyes:
#Overlay-copy-mask w same dims
#swtich to copy - edit mode-
#create mask...

#mrconvert -stride 1,2,3,4 index.mif ../index_strides.nii
#fslroi index_strides.nii index-firstframe_strides.nii 0 1
#fslmaths index-firstframe_strides.nii -thr 1.5 -mas /data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/roi-strides.nii.gz sim/roi_crossings
#fslmaths sim/roi_crossings.nii.gz -bin roi_crossings-bin
#fslmaths index-firstframe_strides.nii  -thr 0.5 -uthr 1.5 -mas roi-strides.nii.gz sim/roi_singles.nii
#fslmaths sim/roi_singles.nii.gz -bin sim/roi_singles_bin.nii.gz

#this is for synthetic inputs
# run in /data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/sim_equation/
# i think i have to make my own masks...
#mrconvert -stride 1,2,3 sim_equation/roi-singles.mif sim_equation/roi_singles.nii
#mrconvert -stride 1,2,3 sim_equation/roi-crossings.mif sim_equation/roi_crossings.nii
# not sure how to save, just a volume for each iteration and merge at the end?

#remmeber to set simulate and sim_input_MTs to True

parser = argparse.ArgumentParser(description='')

parser.add_argument('-o', dest='myresult_output_filename', help='output filename (.npy file). Optional', required=False, nargs=1)
parser.add_argument('-load', dest='myresult_input_filename', help='input filename of previously computed result, to plot it. Optional', required=False, nargs=1)

myargs = parser.parse_args()

iters=100

if (myargs.myresult_input_filename):
    myresult = np.load(myargs.myresult_input_filename)
    print myresult

else:

    for i in range (0,iters):

        #mtr1 = np.zeros([11, 25, 1, iters]) #11 AFD combos, 25 MTR combos, # iterations
        #mtr2 = np.zeros([11, 25, 1, iters])  # 11 AFD combos, 25 MTR combos, # iterations

        #run the sim with 3% noise:
        #mtr_output_iter = '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/sim_equation/sim-equation-iters/mtr_00{}'.format(str(i))

        #subprocess.call(['/data/mril/mril5/ilana/bin/fibermyelin/scripts/fiberMTmyelin_pipeline.py',
                                # #'-mtr', '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/sim_equation/sim-equation-iters/mtr.nii',
                                # '-mtr',mtr_output_iter,
                                # '-Dpar', '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/sim_equation/sim-equation-iters/Dpar.nii',
                                # '-mask', '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/sim_equation/mask.nii.gz',
                                # '-vic', '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/sim_equation/new-icvf.nii',
                                # '-afd', '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/sim_equation/new-afd.nii',
                                # '-afdthresh', '0.1',
                                # '-dirs', '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/sim_equation/new-directions.nii',
                                # '-MTdiff', '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/eddy_corrected_data-on-off-strides.nii',
                                # '-mtr_in','/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/sim_equation/new-mtr.nii',
                                # '-mtw', '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/mt.txt',
                                # '-bvals', '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/nii/on-off.bval',
                                # '-bvecs', '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/nii/eddy_corrected_data-on-off.bvecs',
                                # '-fixel', '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/sim_equation/sim-equation-iters/'])

        #mtr_output_iter = '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/sim_equation/sim-equation-iters-1pc/mtr_00{}'.format(
        #str(i))

        #subprocess.call(['/data/mril/mril5/ilana/bin/fibermyelin/scripts/fiberMTmyelin_pipeline.py',
        #             # '-mtr', '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/sim_equation/sim-equation-iters/mtr.nii',
        #             '-mtr', mtr_output_iter,
        #             '-Dpar',
        #              '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/sim_equation/sim-equation-iters-1pc/Dpar.nii',
        #              '-mask', '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/sim_equation/mask.nii.gz',
        #              '-vic', '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/sim_equation/new-icvf.nii',
        #              '-afd', '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/sim_equation/new-afd.nii',
        #              '-afdthresh', '0.1',
        #              '-dirs',
        #              '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/sim_equation/new-directions.nii',
        #              '-MTdiff',
        #              '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/eddy_corrected_data-on-off-strides.nii',
        #              '-mtr_in', '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/sim_equation/new-mtr.nii',
        #              '-mtw', '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/mt.txt',
        #              '-bvals', '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/nii/on-off.bval',
        #              '-bvecs',
        #              '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/nii/eddy_corrected_data-on-off.bvecs',
        #              '-fixel',
        #              '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/sim_equation/sim-equation-iters-1pc/'])


        mtr_output_iter = '/data_/tardiflab/ilana/mt-diff/qmt_diff_fat_20201201_104549/sim-linear/2pc-FAp7/mtr_00{}'.format(
        str(i))

        subprocess.call(['/data/mril/mril5/ilana/bin/fibermyelin/scripts/fiberMTmyelin_pipeline.py',
                     # '-mtr', '/data_/tardiflab/ilana/mt-diff/mt_diff_64_20191025_104945/sim_equation/sim-equation-iters/mtr.nii',
                     '-mtr', mtr_output_iter,
                     '-Dpar',
                     '/data_/tardiflab/ilana/mt-diff/qmt_diff_fat_20201201_104549/sim-linear/2pc-FAp7/Dpar.nii',
                     '-mask', '/data_/tardiflab/ilana/mt-diff/qmt_diff_fat_20201201_104549/sim/mask.nii.gz',
                     '-afd', '/data_/tardiflab/ilana/mt-diff/qmt_diff_fat_20201201_104549/sim/new-afd.nii',
                     '-afdthresh', '0.1',
                     '-dirs',
                     '/data_/tardiflab/ilana/mt-diff/qmt_diff_fat_20201201_104549/sim/new-directions.nii',
                     '-MTdiff',
                     '/data_/tardiflab/ilana/mt-diff/qmt_diff_fat_20201201_104549/mt-diff-strides.nii',
                     '-mtw', '/data_/tardiflab/ilana/mt-diff/qmt_diff_fat_20201201_104549/mt.txt',
                     '-bvals', '/data_/tardiflab/ilana/mt-diff/qmt_diff_fat_20201201_104549/bvals',
                     '-bvecs','/data_/tardiflab/ilana/mt-diff/qmt_diff_fat_20201201_104549/bvec',
                     '-AD','0.00149',
                     '-RD','0.00038',
                     #'-mtrB0','/data_/tardiflab/ilana/mt-diff/qmt_diff_fat_20201201_104549/sim/mtr-b0.nii',
                     '-fixel',
                     '/data_/tardiflab/ilana/mt-diff/qmt_diff_fat_20201201_104549/sim-linear/2pc-FAp7'])
#break into singles and crossings:
    #subprocess.call(['fslmaths', 'sim-equation-1000/mtr.nii', '-mul', 'roi_singles.nii', 'sim-equation-1000/mtr_singles.nii'])
    #subprocess.call(['fslmaths', 'sim-equation-1000/mtr.nii','-mul', 'roi_crossings.nii', 'sim-equation-1000/mtr_crossings.nii'])



    #get mean and std from fslstats,for singles and crossings:
    #myresult[gap_index, shuffle, 0, 0] = subprocess.check_output(['fslstats', 't1_test_singles.nii.gz', '-M'])
    #myresult[gap_index, shuffle, 0, 1] = subprocess.check_output(['fslstats', 't1_test_crossings.nii.gz', '-M'])
    #myresult[gap_index, shuffle, 1, 0] = subprocess.check_output(['fslstats', 't1_test_singles.nii.gz', '-S'])
    #myresult[gap_index, shuffle, 1, 1] = subprocess.check_output(['fslstats', 't1_test_crossings.nii.gz', '-S'])


    #np.save("myargs.myresult_output_filename", np.array(myresult))

