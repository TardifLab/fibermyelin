#!/usr/bin/python
#/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tues Jun 26 2019
Modified on June 10th  2020 for Inputing Tensors
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

# this is for synthetic inputs
# run in /data_/tardiflab/ilana/midiff_mgh/ir_diff_2p5_20190926_123841/avg-tensor/sim_equation/
# set this in FiberT1Solver.py sim_input_tensor
# just make a mask that covers everything
# not sure how to save, just a volume for each iteration and merge at the end?

parser = argparse.ArgumentParser(description='')

parser.add_argument('-o', dest='myresult_output_filename', help='output filename (.npy file). Optional', required=False,
                    nargs=1)
parser.add_argument('-load', dest='myresult_input_filename',
                    help='input filename of previously computed result, to plot it. Optional', required=False, nargs=1)

myargs = parser.parse_args()

iters = 100

if (myargs.myresult_input_filename):
    myresult = np.load(myargs.myresult_input_filename)
    print myresult

else:

    for i in range(0, iters):

        # run the sim with 3% noise
        tensor_output_iter = '/data_/tardiflab/ilana/midiff_mgh/ir_diff_2p5_20190926_123841/avg-tensor/sim_equation-diff-T1-same-MD/2pc-0pcT1-p8p2-So500/T1_00{}'.format(
            str(i))

        subprocess.call(['/data/mril/mril5/ilana/bin/fibermyelin/scripts/fibermyelin_pipeline.py',
                         '-t1', tensor_output_iter,
                         '-Dpar',
                         '/data_/tardiflab/ilana/midiff_mgh/ir_diff_2p5_20190926_123841/avg-tensor/sim_equation-diff-T1-same-MD/2pc-0pcT1-p8p2-So500/Dpar.nii',
                         '-mask', '/data_/tardiflab/ilana/midiff_mgh/ir_diff_2p5_20190926_123841/avg-tensor/sim_equation-diff-T1/mask.nii',
                         '-vic', '/data_/tardiflab/ilana/midiff_mgh/ir_diff_2p5_20190926_123841/hardi-analysis/AMICO/NODDI/FIT_ICVF-strides.nii.gz',
                         '-afd', '/data_/tardiflab/ilana/midiff_mgh/ir_diff_2p5_20190926_123841/avg-tensor/sim_equation-diff-T1/new-afd-p8p2.nii',
                         '-afdthresh', '0.1',
                         '-dirs',
                         '/data_/tardiflab/ilana/midiff_mgh/ir_diff_2p5_20190926_123841/avg-tensor/sim_equation-diff-T1/new-directions.nii',
                         '-IRdiff',
                         '/data_/tardiflab/ilana/midiff_mgh/ir_diff_2p5_20190926_123841/ir-analysis/ir-diff-strides.nii',
                         '-TIs', '/data_/tardiflab/ilana/midiff_mgh/ir_diff_2p5_20190926_123841/TIcomb.txt',
                         '-ADin',
                         '/data_/tardiflab/ilana/midiff_mgh/ir_diff_2p5_20190926_123841/avg-tensor/sim_equation-diff-T1-same-MD/new-AD.nii',
                         '-RDin',
                         '/data_/tardiflab/ilana/midiff_mgh/ir_diff_2p5_20190926_123841/avg-tensor/sim_equation-diff-T1-same-MD/new-RD.nii',
                         '-bvals', '/data_/tardiflab/ilana/midiff_mgh/30bvals-4avg-rmb0.txt',
                         '-bvecs',
                         '/data_/tardiflab/ilana/midiff_mgh/30bvecs-4avg-rmb0.txt',
                         '-TR', '3000',
                         '-TE', '56',
                         '-fixel',
                         '/data_/tardiflab/ilana/midiff_mgh/ir_diff_2p5_20190926_123841/avg-tensor/sim_equation-diff-T1-same-MD/2pc-0pcT1-p8p2-So500/'])

