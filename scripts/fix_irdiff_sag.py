#!/usr/bin/python
#/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  12 12:12 2018

@author: jcampbel
"""

import numpy as np
import nibabel as nib
import argparse

IR_diff_img=nib.load("unshuffle_diff_series4_tmp.nii")

#y is stored in first voxel index and we need to flip that:

new_array=np.zeros(IR_diff_img.header.get_data_shape())

for y in range(IR_diff_img.header.get_data_shape()[0]):
    for z in range(IR_diff_img.header.get_data_shape()[1]):
        for x in range(IR_diff_img.header.get_data_shape()[2]):
            for t in range(IR_diff_img.header.get_data_shape()[3]):
                #new_array[y,z,x,t]=IR_diff_img.get_data()[IR_diff_img.header.get_data_shape()[0]-1-y,z,x,t]
                #test: flip x:
                new_array[y,z,x,t]=IR_diff_img.get_data()[y,z,IR_diff_img.header.get_data_shape()[2]-1-x,t]                
                
new_IR_diff_img = nib.Nifti1Image(new_array, IR_diff_img.affine, IR_diff_img.header)
nib.save(new_IR_diff_img, "unshuffle_diff_series4_tmp_xflip.nii")