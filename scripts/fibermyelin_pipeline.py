#!/usr/bin/python
#/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 13:43:18 2017

@author: jcampbel
"""

'''This is a pipeline to process IR-diffusion data'''


#from pyminc.volumes.factory import *
from fiberT1.FiberT1Solver import FiberT1Solver
import os
import numpy as np
import nibabel as nib
from dipy.io import read_bvals_bvecs
from dipy.core.gradients import gradient_table
import vtk
#import sys
import argparse
import os.path
from os import path

#constants:
EPSILON=9E-9

#hardcoded stuff:

#comparison for -sort option:
xz=False
xy=True



#file orientation (if False, assumes axial)
global sagittal
sagittal=False

global axial

if (sagittal==False):
    axial=True

global set_Dpar_equal
set_Dpar_equal=True #HERE this has to be true right now have to set in FiberT1Solver.py too

global bedpostx
bedpostx=False

global just_b0
just_b0=False #have to set in FiberT1Solver.py too

# Voxel T1 averaging settings
global true_afd_thresh  # used to identify sub-afdthresh fibers that will have their T1 times replaced after fitting 
if just_b0:
    true_afd_thresh = 0.0  # needs to be zero so that T1 is calculated in each voxel 
    print('')
    print('***********************')
    print("true_afd_thresh = %0.2f" % true_afd_thresh)
    print('***********************')
    print('')
else:
    true_afd_thresh = 0.2
    print('')
    print('***********************')
    print("true_afd_thresh = %0.2f" % true_afd_thresh)
    print('***********************')
    print('')

global on_true_afd_thresh # activates within-voxel fiber T1 averaging 
on_true_afd_thresh = False
# NOTE - when running fibermyelin_pipeline.py, set afdthresh = 0.0 so that T1 is calculated for ALL fibers in a voxel
# fibermyelin_pipeline treats within-voxel T1 times in 1 of 2 ways (depending on fiber AFD):
# 1) if any fibers in the voxel are super-true_afd_thresh, T1 times of sub-true_afd_tresh fibers are set to the weighted average of super-thresh fiber T1 times
# 2) if only sub-true_afd_thresh fibers in the voxel, these fiber T1 times are calculated using the SE-IR equation 
# Outcome - same final T1 image/fixel dimensions as initial afd/directions image/fixel dimensions, smoother tract T1 maps

#for visualization W/L:
#asparagus 500-800,1500
#human brain 550-575-600,725-750-800
global vis_min
vis_min=500
global vis_max
vis_max=1000
global vis_range
vis_range=vis_max-vis_min


#end hardcoded stuff

if just_b0:
    set_Dpar_equal=False


print('This a script to process IR-diffusion data\n\
NB: it does not currently handle data acquired with angulation\n\
NB: it assumes the following has already been done:\n\
    -registration of all images and sampling to same voxel space\n\
    -AFD computation with MRtrix\n\
    -optionally, tensor and/or NODDI fit\n\
NB: there is an option to hardcode vic in FiberT1Solver; you must input a vic file regardless:\n\
    just input your mask. Edit FiberT1Solver global var fix_vic to actually use vic.\n\
NB: using vic assumes the same tortuosity for all fibers; we are thinking about an alternative.\n\
NB: there is currently a Johnson noise term added to all images.\n\
NB: the data are assumed to have been acquired with diffusion fastest varying, and unshuffled to TI fastest varying.\n\
    i.e., bvals/bvecs have diffusion fastest varying, and input IRdiff images do not.\n\
NB: the visualization W/L is a WIP\n\
NB: initial parameter values are hardcoded in FiberT1Solver::GetT1s,\n\
    as is the fitting algorithm: currently trust region reflective (\'trf\').\n\
    Upper and lower bounds are also hardcoded\n\
    Comment/uncomment appropriate lines to change fitting algorithm and fitting options.\n\
NB: the b=0 images are currently weighted significantly more than the rest.\n\
NB: future additions that aren\'t here: CSF compartment, GM compartment\n\
NB: the fixel directory output is mandatory and only works properly for axial images.')

#DO: add fixel directories or other output format for thresholded AFDs (and possibly Dpar)

parser = argparse.ArgumentParser(description='')

parser.add_argument('-mask', dest='mask_image_filename', help='input mask filename, for computation or visualization. If for visualization, must be the same mask previously used for computation', required=True, nargs=1)
parser.add_argument('-vic', dest='vic_image_filename', help='input vic filename, for computation only (not neede for mean tensor option).', required=False, nargs=1)
parser.add_argument('-S0', dest='S0_image_filename', help='input S0 filename, optinally used to stabilize fits', required=False, nargs=1)
parser.add_argument('-T1dw', dest='T1dw_image_filename', help='input computed average diffusion-weighted T1, optinally used as initial guess in fits', required=False, nargs=1)
parser.add_argument('-just_T1dw', dest='just_T1dw', help='flag to compute only diffusion-weighted T1 (direction averaged)', required=False, action='store_true')
parser.add_argument('-afd', dest='afd_image_filename', help='input AFD filename, currently required even for -vis', required=True, nargs=1)
parser.add_argument('-afdthresh', dest='AFD_thresh', help='input AFD threshold, required.  For -vis option, AFD threshold must be that used to compute the pre-computed T1 map', required=True, nargs=1)
parser.add_argument('-dirs', dest='dirs_image_filename', help='input directions filename, currently required even for -vis', required=True, nargs=1)
parser.add_argument('-IRdiff', dest='IR_diff_image_filename', help='input IR diffusion weighted image series; required for full computation',  required=False, nargs=1)
parser.add_argument('-TIs', dest='TIs_filename', help='input text file of TIs for each slice; required for full computation',  required=False, nargs=1)
parser.add_argument('-bvals', dest='bvals_filename', help='input text file of bvals; required for full computation',  required=False, nargs=1)
parser.add_argument('-bvecs', dest='bvecs_filename', help='input text file of bvecs; required for full computation',  required=False, nargs=1)
parser.add_argument('-TR', dest='TR', help='input nominal TR for short-TR steady state equation (ms)',required=True,nargs=1)
parser.add_argument('-TE', dest='TE', help='input nominal TE for short-TR steady state equation (ms)',required=True,nargs=1)
parser.add_argument('-viso', dest='viso_image_filename', help='input viso filename, for computation only. Optional', required=False, nargs=1)
parser.add_argument('-B1', dest='B1_image_filename', help='input B1 filename, for computation only. Optional', required=False, nargs=1)
parser.add_argument('-t1', dest='t1_image_filename', help='T1 filename: output filename if doing full computation; pre-computed T1 filename if -vis option is selected', required=True, nargs=1)
parser.add_argument('-S0_out', dest='S0out_filename', help='S0-out filename: output S0 image filename', required=False, nargs=1)
parser.add_argument('-cost', dest='cost_filename', help='cost filename: output cost image filename, to see quality of fit', required=False, nargs=1)
parser.add_argument('-fixel',dest='fixel_dir_name', help='output T1 fixel directory name; currently implemented only for axial images!!', required=True, nargs=1)
parser.add_argument('-Dpar', dest='Dpar_image_filename', help='output D_parallel filename', required=False, nargs=1)
parser.add_argument('-ADin', dest='AD_image_filename', help='Optional input AD (axial diffusivity) for simulations, make sure sim_input_tensor is set', required=False, nargs=1)
parser.add_argument('-RDin', dest='RD_image_filename', help='Optional input RD (radial diffusivity) for simulations, make sure sim_input_tensor is set, need to use with ADin', required=False, nargs=1)
parser.add_argument('-vis', dest='visualize', help='visualize previously computed fiber T1s', required=False, action='store_true')
parser.add_argument('-sort', dest='sortT1', help='print out previously computed fiber T1s in 2-fiber voxels, sorted by orientation', required=False, action='store_true')

myargs = parser.parse_args()


#currently reads in .nii files
#assumes already done:
#standard registration? be careful because different slices have different TIs, etc.
#everything has to be sampled the same way, and like the IR-diff data, i.e., that number of slices. don't change the slicing because it determines TI
#i.e. check whether the AFD, etc. have the same sampling and fix them if not.

#steady-state (non-infinite TR) equation is always used
#it requires TR and TE to be input
#however, it won't work if slices shuffled and TR short enough to matter; we are handling this with dummies.

#I believe the bvecs (from dcm2nii) are in voxel space with strides -1,2,3. That is how the code is right now.

#NOTE: the bvecs and bvals have to be edited to be correct with the dummy scan removed.

#steps: afd, directions
#convert the AFD and directions files to non-sparse format:
#fixel2voxel <afd.mif> split_data <afd_voxel.mif> ...
#for our pipeline: fixel2voxel fixel_dir/afd.mif split_data afd_voxel.mif
#fixel2voxel <directions.mif> split_dir <directions_voxel.mif> ...
#for our pipeline: fixel2voxel fixel_dir/directions.mif split_dir directions_voxel.mif
#-crop if necessary, mrcrop...0 offset

#for axial:#mrconvert -stride 1,2,3,4 <afd> ...
#mrconvert -stride 1,2,3,4 afd_voxel.mif afd_voxel_strides.nii
#mrconvert -stride 1,2,3,4 directions_voxel.mif directions_voxel_strides.nii


#for sagittal:
#mrconvert -stride 3,1,2 <file.mif> <file.nii>

#steps: ir-diff
#make sure it came directly from dcm
#axial:
#need to mrconvert -stride 1,2,3,4
#
#if sagittal, mrconvert -stride 3,1,2,4
#then, **flip in x, using fix_IRdiff_sag.py** OR just set strides to -3,1,2,4
#I think this is because of something in Ilana's unshuffling script, and we could change that and not flip here

#do same to mask as to files that match it (first 3 dims)
#e.g., for axial, I needed to set strides to 1,2,3
#for sagittal, 3,1,2

#NB!! NODDI processing erases the strides.  So, what I did was resample the file, reversing the negative stride direction (x)
#it might work to use strides -1,2,3, but I haven't tested it.
#I used script fix_NODDI_trans

if (myargs.visualize or myargs.sortT1):
    T1_array=nib.load(myargs.t1_image_filename[0]).get_data()
else:
    if path.exists(myargs.fixel_dir_name[0]):
        print ("output directory %s exists, overwriting content\n", myargs.fixel_dir_name[0])
    else:
        os.mkdir(myargs.fixel_dir_name[0])

    IR_diff_img = nib.Nifti2Image.from_image(nib.load(myargs.IR_diff_image_filename[0]))
    if myargs.vic_image_filename:
        vic_img = nib.Nifti2Image.from_image(nib.load(myargs.vic_image_filename[0]))
    if myargs.viso_image_filename:
        viso_img=nib.Nifti2Image.from_image(nib.load(myargs.viso_image_filename[0]))
    if myargs.S0_image_filename:
        S0_img=nib.Nifti2Image.from_image(nib.load(myargs.S0_image_filename[0]))
    if myargs.T1dw_image_filename:
        T1dw_img=nib.Nifti2Image.from_image(nib.load(myargs.T1dw_image_filename[0]))
    if myargs.B1_image_filename:
        B1_img=nib.Nifti2Image.from_image(nib.load(myargs.B1_image_filename[0]))

    if myargs.AD_image_filename:
        AD_img=nib.Nifti2Image.from_image(nib.load(myargs.AD_image_filename[0]))
    if myargs.RD_image_filename:
        RD_img=nib.Nifti2Image.from_image(nib.load(myargs.RD_image_filename[0]))

#For visualize, sort, and full computation:  
#AFD_img=nib.load(myargs.afd_image_filename[0])    
AFD_img=nib.Nifti2Image.from_image(nib.load(myargs.afd_image_filename[0]))    
max_fibers=AFD_img.header.get_data_shape()[3]    



#fiber_dirs_img=nib.load(myargs.dirs_image_filename[0]) 
fiber_dirs_img=nib.Nifti2Image.from_image(nib.load(myargs.dirs_image_filename[0])) 



#FSL BEDPOSTX version:

if (bedpostx):
    max_fibers=2
    AFD_img_fsl=[]
    fiber_dirs_img_fsl=[]
    AFD_img_fsl.append(nib.load("bedpost_outputs_strides/mean_f1samples_strides.nii"))
    AFD_img_fsl.append(nib.load("bedpost_outputs_strides/mean_f2samples_strides.nii"))
    fiber_dirs_img_fsl.append(nib.load("bedpost_outputs_strides/dyads1_strides.nii"))
    fiber_dirs_img_fsl.append(nib.load("bedpost_outputs_strides/dyads2_strides.nii"))


#mask_img=nib.load(myargs.mask_image_filename[0])    
mask_img=nib.Nifti2Image.from_image(nib.load(myargs.mask_image_filename[0]))    
    
voxels=np.where(mask_img.get_data()>=EPSILON)  

number_of_fibers_array=np.zeros(mask_img.header.get_data_shape(),int)

#make a new AFD and fiber_dirs array for thresholded fiber dirs (they aren't sequential,at least beyond first)
AFD_array=np.zeros(AFD_img.header.get_data_shape())
fiber_dirs_array=np.zeros(fiber_dirs_img.header.get_data_shape())

#FSL BEDPOSTX version:
if (bedpostx):

    AFD_array=np.zeros([AFD_img_fsl[0].header.get_data_shape()[0], AFD_img_fsl[0].header.get_data_shape()[1], AFD_img_fsl[0].header.get_data_shape()[2], max_fibers])
    fiber_dirs_array=np.zeros([AFD_img_fsl[0].header.get_data_shape()[0], AFD_img_fsl[0].header.get_data_shape()[1], AFD_img_fsl[0].header.get_data_shape()[2], 3*max_fibers])


AFD_thresh=float(myargs.AFD_thresh[0])

TR=float(myargs.TR[0])
TE=float(myargs.TE[0])




voxelcounter=0

for j in range(len(voxels[0])):

    #get number of super-threshold AFDs=number_of_fibers
    number_of_fibers=0
    for l in range(max_fibers):
        #ASSUMING ORDERED and stopping at 3:
        if (l<3):
            if AFD_img.get_data()[voxels[0][j],voxels[1][j],voxels[2][j],l]>AFD_thresh:
                AFD_array[voxels[0][j],voxels[1][j],voxels[2][j],number_of_fibers]=AFD_img.get_data()[voxels[0][j],voxels[1][j],voxels[2][j],l]
                fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],3*number_of_fibers]=fiber_dirs_img.get_data()[voxels[0][j],voxels[1][j],voxels[2][j],3*l]
                fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],3*number_of_fibers+1]=fiber_dirs_img.get_data()[voxels[0][j],voxels[1][j],voxels[2][j],3*l+1]
                fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],3*number_of_fibers+2]=fiber_dirs_img.get_data()[voxels[0][j],voxels[1][j],voxels[2][j],3*l+2]
                number_of_fibers+=1


    number_of_fibers_array[voxels[0][j],voxels[1][j],voxels[2][j]]=number_of_fibers

#FSL BEDPOSTX version:
#note that fiber dirs in this case are in "strided" voxel space
if (bedpostx):
    for j in range(len(voxels[0])):
        for l in range(max_fibers):
            AFD_array[voxels[0][j],voxels[1][j],voxels[2][j],l]=AFD_img_fsl[l].get_data()[voxels[0][j],voxels[1][j],voxels[2][j]]
            fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],3*l]=-1.0*fiber_dirs_img_fsl[l].get_data()[voxels[0][j],voxels[1][j],voxels[2][j],0]
            fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],3*l+1]=fiber_dirs_img_fsl[l].get_data()[voxels[0][j],voxels[1][j],voxels[2][j],1]
            fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],3*l+2]=fiber_dirs_img_fsl[l].get_data()[voxels[0][j],voxels[1][j],voxels[2][j],2]
        number_of_fibers_array[voxels[0][j],voxels[1][j],voxels[2][j]]=max_fibers

if not (myargs.visualize or myargs.sortT1):


    number_of_slices=0
    with open(myargs.TIs_filename[0],'r') as f:
        for line in f:
            number_of_TIs = len(line.split())
            number_of_slices+=1

    TIs_allslices=np.zeros([number_of_slices,number_of_TIs])
    #TRs_allslices=np.zeros([number_of_slices,number_of_TIs])


    i=0
    with open(myargs.TIs_filename[0],'r') as f:
        for line in f:
            TIs_allslices[i]=np.array(line.split()).astype(np.float)
            i+=1




    #get the bvals and bvecs for this acquisition (currently assuming same for each TI!)
    #test2: bvals, bvecs = read_bvals_bvecs("bvals", "bvecs")
    bvals, bvecs = read_bvals_bvecs(myargs.bvals_filename[0],myargs.bvecs_filename[0])
    grad_table = gradient_table(bvals, bvecs)

    #the bvecs are in strided voxel space, so we negate here  x (2) and y (0)
    if (sagittal): #DO: check this again
        grad_table.bvecs[:,0]=-1.0*grad_table.bvecs[:,0]
        grad_table.bvecs[:,2]=-1.0*grad_table.bvecs[:,2]
    else: #axial (even if prone)
        grad_table.bvecs[:,0]=-1.0*grad_table.bvecs[:,0]
    #print grad_table.bvecs[:,0]



    #number of DWIs (assume one b=0, at beginning). assume x,y,z,diff

    number_of_DWIs=IR_diff_img.get_data().shape[3]



    T1_array=np.zeros(AFD_img.header.get_data_shape())
    Dparfit_array=np.zeros(AFD_img.header.get_data_shape())

    if myargs.S0out_filename:
        S0out_array = np.zeros(mask_img.header.get_data_shape())
    if myargs.cost_filename:
        cost_array = np.zeros(mask_img.header.get_data_shape())

    for j in range(len(voxels[0])):

        voxelcounter += 1

        number_of_fibers=number_of_fibers_array[voxels[0][j],voxels[1][j],voxels[2][j]]

        AFD=np.zeros(number_of_fibers,float)
        fiber_dirs=np.zeros((number_of_fibers,3),float)


        for i in range(0,number_of_fibers):
            #for unthresholded only:
            #AFD[i]=AFD_img.get_data()[voxels[0][j],voxels[1][j],voxels[2][j],i]
#            fiber_dirs[i,0]=fiber_dirs_img.get_data()[voxels[0][j],voxels[1][j],voxels[2][j],3*i]
#            fiber_dirs[i,1]=fiber_dirs_img.get_data()[voxels[0][j],voxels[1][j],voxels[2][j],3*i+1]
#            fiber_dirs[i,2]=fiber_dirs_img.get_data()[voxels[0][j],voxels[1][j],voxels[2][j],3*i+2]
            AFD[i]=AFD_array[voxels[0][j],voxels[1][j],voxels[2][j],i]


            fiber_dirs[i,0]=fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],3*i]
            fiber_dirs[i,1]=fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],3*i+1]
            fiber_dirs[i,2]=fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],3*i+2]

        print('\n---Fit [%i/%i]' % (voxelcounter, len(voxels[0])))
        print('number_of_fibers %i:' % number_of_fibers)
        print(fiber_dirs)
        print("AFDs:")
        print(AFD)

        if myargs.vic_image_filename:
            vic=vic_img.get_data()[voxels[0][j],voxels[1][j],voxels[2][j]]
            print('vic: %f' % viso)
        else:
            vic=0

        if myargs.viso_image_filename:
            viso=viso_img.get_data()[voxels[0][j],voxels[1][j],voxels[2][j]]
            print('viso: %f' % viso)
        else:
            viso=0

        if myargs.B1_image_filename:
            B1=B1_img.get_data()[voxels[0][j],voxels[1][j],voxels[2][j]]
        else:
            B1=1.0
        if myargs.S0_image_filename:
            S0=S0_img.get_data()[voxels[0][j],voxels[1][j],voxels[2][j]]
            print('S0: %f' % S0)
        else:
            S0=500
        if myargs.T1dw_image_filename:
            T1dw=T1dw_img.get_data()[voxels[0][j],voxels[1][j],voxels[2][j],0]
            print('Avg diff-weighted T1: %f' % T1dw)
        else:
            T1dw=750 #the initial value for T1

        if myargs.AD_image_filename: #assume both are provided
            AD = AD_img.get_data()[voxels[0][j], voxels[1][j], voxels[2][j],:]
            RD = RD_img.get_data()[voxels[0][j], voxels[1][j], voxels[2][j],:]
        else:
            AD = 0
            RD = 0

        IR_DWIs=np.zeros((number_of_TIs, number_of_DWIs),float)

        for i in range(0,number_of_TIs):
            IR_DWIs=IR_diff_img.get_data()[voxels[0][j],voxels[1][j],voxels[2][j],:]

        #TIs is the TIs for this slice coordinate  **if axial, slices in z**
        #**if sagittal, slices in x, which is still stored in dim 2**


        TIs=TIs_allslices[voxels[2][j]]

        if myargs.just_T1dw:
            set_Dpar_equal = False


        #TRs=TRs_allslices[voxels[2][j]] #HERE do this when have the whole TR history in signal comp.

        S0out = 0
        cost = 0
        if (number_of_fibers>0):
            t1solver=FiberT1Solver()

            #to test fixing D, we have to giev it the voxel coord

            t1solver.SetInputData(fiber_dirs,AFD,IR_DWIs,TIs,grad_table,vic,S0,T1dw,TR,TE,voxels[0][j],voxels[1][j],voxels[2][j],sagittal,set_Dpar_equal,viso,B1,AD,RD,myargs.just_T1dw)
            T1s=np.zeros(number_of_fibers)
            Dparfit=np.zeros(number_of_fibers)


            T1sandDparfit=np.zeros([number_of_fibers,2])
            FullFit =t1solver.GetT1s()
            if (FullFit is not None): #fit success
                T1sandDparfit = FullFit.x
                #T1sandDparfit=t1solver.GetT1s()
                if myargs.S0out_filename:
                    if (just_b0 or myargs.just_T1dw):
                        S0out = T1sandDparfit[2]
                    else:
                        S0out = T1sandDparfit[number_of_fibers+2]
                if myargs.cost_filename:
                    cost = FullFit.cost
            else:
                T1sandDparfit = None


        #hack for just_b0 just_T1dw:
        if (just_b0 or myargs.just_T1dw):
            if(number_of_fibers==0):
                number_of_fibers = 1
                T1s=np.zeros(number_of_fibers)
                Dparfit=np.zeros(number_of_fibers)
                T1sandDparfit=np.zeros(number_of_fibers+1)

        for i in range(0,number_of_fibers):
            if (T1sandDparfit is not None):
                if (set_Dpar_equal):

                    T1s[i]=T1sandDparfit[i]
                    Dparfit[i]=T1sandDparfit[number_of_fibers]
                    ## IRL commented this whole section because I'm not sure what's going on...
                    #                     if on_true_afd_thresh: # do the voxel T1 averaging for sub afdthresh fibers in voxels with >=1 superthresh fibers, fiber T1s in voxels with only subthresh fibers are just_b0 T1s already
                    #
                    #                         #if len(AFD)==1: # Single fiber voxel  //   not needed since the just_b0 equation is used in voxels with only subthresh fibers
                    # 		        #    if AFD[i]<true_afd_thresh:
                    # 			#        T1s[0] = 0.0 # assigning nan doesn't work currently. Set to 0 for now since this is a minority case, exclude this voxel manually (as of Dec. 4, 2019)
                    #
                    #                         if len(AFD)==2: # Two fibers
                    # 		            if (AFD[0]<true_afd_thresh and AFD[1]>true_afd_thresh):   # LR
                    # 			        T1s[0] = T1sandDparfit[1]
                    # 		            elif (AFD[0]>true_afd_thresh and AFD[1]<true_afd_thresh): # RL
                    # 			        T1s[1] = T1sandDparfit[0]
                    # #Added Dec. 10, 2019
                    # 		            elif AFD[0]<true_afd_thresh and AFD[1]<true_afd_thresh: # LL
                    #                                 T1s[0] = T1sandDparfit[1] # set the T1 of the first to the same as the last fiber
                    #
                    # # Original before Dec .10, 2019
                    # 		            #elif AFD[0]<true_afd_thresh and AFD[1]<true_afd_thresh: # LLL
                    #                             #    mean_voxel_T1 = (AFD[0]*T1s[0] + AFD[1]*T1s[1])/(AFD[0] + AFD[1])
                    # 			    #	T1s[0] = mean_voxel_T1 # Set the T1 of each fiber to the weighted average T1 of all the fibers
                    # 			    #   T1s[1] = mean_voxel_T1
                    #
                    #                         if len(AFD)==3: # Three fibers (current maximum n_fibers per voxel is 3 - Dec. 3, 2019)
                    # 		            if (AFD[0]<true_afd_thresh and AFD[1]<true_afd_thresh and AFD[2]>true_afd_thresh):   # LLR ; L == sub-afdthresh fiber, R == super-afdthresh fiber. Order of fibers in array must be considered.
                    # 			        T1s[0] = T1sandDparfit[2]
                    # 			        T1s[1] = T1sandDparfit[2]
                    # 		            elif (AFD[0]<true_afd_thresh and AFD[1]>true_afd_thresh and AFD[2]>true_afd_thresh): # LRR
                    # 			        T1s[0] = (AFD[1]*T1sandDparfit[1] + AFD[2]*T1sandDparfit[2])/(AFD[1]+AFD[2])
                    # 		            elif (AFD[0]>true_afd_thresh and AFD[1]>true_afd_thresh and AFD[2]<true_afd_thresh): # RRL
                    # 			        T1s[2] = (AFD[0]*T1sandDparfit[0] + AFD[1]*T1sandDparfit[1])/(AFD[0]+AFD[1])
                    # 		            elif (AFD[0]>true_afd_thresh and AFD[1]<true_afd_thresh and AFD[2]<true_afd_thresh): # RLL
                    # 			        T1s[1] = T1sandDparfit[0]
                    # 			        T1s[2] = T1sandDparfit[0]
                    # 		            elif (AFD[0]<true_afd_thresh and AFD[1]>true_afd_thresh and AFD[2]<true_afd_thresh): # LRL
                    # 			        T1s[0] = T1sandDparfit[1]
                    # 			        T1s[2] = T1sandDparfit[1]
                    # 		            elif (AFD[0]>true_afd_thresh and AFD[1]<true_afd_thresh and AFD[2]>true_afd_thresh): # RLR
                    # 			        T1s[1] = (AFD[0]*T1sandDparfit[0] + AFD[2]*T1sandDparfit[2])/(AFD[0]+AFD[2])
                    # #Added Dec. 10, 2019
                    # 		            elif (AFD[0]<true_afd_thresh and AFD[1]<true_afd_thresh and AFD[2]<true_afd_thresh): # LLL
                    #                                  T1s[0] = T1sandDparfit[2] # Set the T1 of each fiber to the T1 of the final counted fiber in the voxel
                    # 			         T1s[1] = T1sandDparfit[2] # Required because T1 is only fit for the final fiber using the SE-IR equation (related to looping structure in FiberT1Solver.py ~line 744)
                    # 			                                   # No diffusion term used, all fibers initiliazed with the same T1 values
                    #
                    #
                    # # Original before Dec. 10, 2019
                    # 		            #elif AFD[0]<true_afd_thresh and AFD[1]<true_afd_thresh and AFD[2]<true_afd_thresh: # LLL
                    #                             #    mean_voxel_T1 = (AFD[0]*T1s[0] + AFD[1]*T1s[1] + AFD[2]*T1s[2])/(AFD[0] + AFD[1] + AFD[2])
                    # 			    #	T1s[0] = mean_voxel_T1 # Set the T1 of each fiber to the weighted average T1 of all the fibers
                    # 			    #   T1s[1] = mean_voxel_T1
                    # 			    #   T1s[2] = mean_voxel_T1
                    ## IRL end of commented section
                    print('Fiber T1: %d' % T1s[i])
                    print('Voxel Dpar: %f' % Dparfit[i])
                elif (just_b0 or myargs.just_T1dw):
                    T1s[i]=T1sandDparfit[0]
                    print('T1: %d' % T1s[0])
                else:
                    T1s[i]=T1sandDparfit[2*i]
                    Dparfit[i]=T1sandDparfit[2*i+1]
                    print('T1: %d' % T1s[i])
                    print('Dpar: %f' % Dparfit[i])

                T1_array[voxels[0][j], voxels[1][j], voxels[2][j], i] = T1s[i]           # Store the T1 times in the array 
                Dparfit_array[voxels[0][j], voxels[1][j], voxels[2][j], i] = Dparfit[i]

                if myargs.S0out_filename:
                    S0out_array[voxels[0][j], voxels[1][j], voxels[2][j]] = S0out
                if myargs.cost_filename:
                    cost_array[voxels[0][j], voxels[1][j], voxels[2][j]] = cost
            else: #fit is none, set T1 to B=0 value

                T1_array[voxels[0][j], voxels[1][j], voxels[2][j], i] = T1dw

        #optional: this plots the single voxel right now:
        #t1solver.VisualizeFibers()


#        if (voxelcounter%1000 == 0):
            #output the T1 map.
            #start as a nonsparse copy of AFD file
    T1_img = nib.Nifti2Image(T1_array, AFD_img.affine, AFD_img.header)
    nib.save(T1_img, myargs.t1_image_filename[0])


    Dparfit_img = nib.Nifti2Image(Dparfit_array, AFD_img.affine, AFD_img.header)
    if myargs.Dpar_image_filename is None:
        nib.save(Dparfit_img, "Dpar.nii")
    else:
        nib.save(Dparfit_img, myargs.Dpar_image_filename[0])

    if myargs.S0out_filename:
        S0out_img = nib.Nifti2Image(S0out_array, mask_img.affine, mask_img.header)
        nib.save(S0out_img, myargs.S0out_filename[0])
    if myargs.cost_filename:
        cost_img = nib.Nifti2Image(cost_array, mask_img.affine, mask_img.header)
        nib.save(cost_img, myargs.cost_filename[0])

    print(' ') 
    print('**************************************************************')
    print('Parameter calculations complete. Saving images and fixel files.')
    print('**************************************************************')

    T1_array_zeroed=np.zeros(AFD_img.header.get_data_shape())

    for j in range(len(voxels[0])):  # Do not assign zero values here - this is done earlier 
        counter=0
        for l in range(max_fibers):
            # Since afdthresh = 0.0, there is always a fiber / should really be an index threshold if done this way
            if AFD_img.get_data()[voxels[0][j],voxels[1][j],voxels[2][j],l]>AFD_thresh:  # if fiber is present at all, count it (since afdthresh is fiber presence threshold, not fit threshold)
                T1_array_zeroed[voxels[0][j],voxels[1][j],voxels[2][j],l]=T1_array[voxels[0][j],voxels[1][j],voxels[2][j],counter]
                counter+=1
# Original
#            if AFD_img.get_data()[voxels[0][j],voxels[1][j],voxels[2][j],l]>AFD_thresh:  # if fiber is present at all, count it (since afdthresh is fiber presence threshold, not fit threshold)
#                T1_array_zeroed[voxels[0][j],voxels[1][j],voxels[2][j],l]=T1_array[voxels[0][j],voxels[1][j],voxels[2][j],counter]
#                counter+=1
#            else:
#                T1_array_zeroed[voxels[0][j],voxels[1][j],voxels[2][j],l]=0.0  # T1 values of subthresh afd fibers are set to zero // ignored since the afdthresh input at runtime is absolute thresh for fiber presence 
#	        print('***********************************************************')
#	        print('Fiber below afdthresh!!! Zero T1 value assigned to fiber!!!')
#	        print('***********************************************************')

    if (False):
        #output the T1s to match the directions initially input:
        #the intent is to then use voxel2fixel, but it fails.  It writes the first T1 to all fixels.
        T1_img_zeroed = nib.Nifti2Image(T1_array_zeroed, AFD_img.affine, AFD_img.header)
        nib.save(T1_img_zeroed, "t1_zeroed.nii") # the zeroed image was never the one being saved - the

    #output the directions, t1, and index files in sparse format.
    #I'm doing this for the thresholded data
    #index is the size of AFD, first frame


    new_size=[AFD_img.header.get_data_shape()[0],AFD_img.header.get_data_shape()[1], AFD_img.header.get_data_shape()[2],2]

    index_array=np.zeros(new_size)

    #temporarily make these as enormous as possible:


    #thresh_dirs_array=np.zeros([AFD_img.header.get_data_shape()[0]*AFD_img.header.get_data_shape()[1]*AFD_img.header.get_data_shape()[2],3,1])
    thresh_dirs_array = np.zeros([len(voxels[0]) * max_fibers, 3, 1])
    #t1_fixel_array=np.zeros([AFD_img.header.get_data_shape()[0]*AFD_img.header.get_data_shape()[1]*AFD_img.header.get_data_shape()[2],1,1])
    t1_fixel_array = np.zeros([len(voxels[0]) * max_fibers, 3, 1])

    counter=0
    #for i in range(AFD_img.header.get_data_shape()[0]*AFD_img.header.get_data_shape()[1]*AFD_img.header.get_data_shape()[2])
    for i in range(len(voxels[0])):#for just the mask
        fibercounter=0
        for l in range(max_fibers):
            #ASSUMING ORDERED and stopping at 3:
            if (l<3):
                if AFD_img.get_data()[voxels[0][i],voxels[1][i],voxels[2][i],l]>AFD_thresh:  # will always run if the fiber exists (afdthresh right now is the fiber existence threshold)
                    if (index_array[voxels[0][i],voxels[1][i],voxels[2][i],0]==0):#set on first fiber

                        index_array[voxels[0][i],voxels[1][i],voxels[2][i],1]=counter
                    thresh_dirs_array[counter,0,0]=fiber_dirs_img.get_data()[voxels[0][i],voxels[1][i],voxels[2][i],3*l]
                    thresh_dirs_array[counter,1,0]=fiber_dirs_img.get_data()[voxels[0][i],voxels[1][i],voxels[2][i],3*l+1]
                    thresh_dirs_array[counter,2,0]=fiber_dirs_img.get_data()[voxels[0][i],voxels[1][i],voxels[2][i],3*l+2]

		    if (just_b0 or myargs.just_T1dw):
		        t1_fixel_array[counter,0,0]=T1_array_zeroed[voxels[0][i],voxels[1][i],voxels[2][i],0]
			#print("%d\n" % t1_fixel_array[counter,0,0])
		    else:
                        t1_fixel_array[counter,0,0]=T1_array_zeroed[voxels[0][i],voxels[1][i],voxels[2][i],l]

                    counter+=1
                    fibercounter+=1
                    index_array[voxels[0][i],voxels[1][i],voxels[2][i],0]=fibercounter


    if not os.path.exists(myargs.fixel_dir_name[0]):
        os.makedirs(myargs.fixel_dir_name[0])

               
    index_img = nib.Nifti2Image(index_array,AFD_img.affine)
    
    
    nib.save(index_img, myargs.fixel_dir_name[0]+"/index.nii")
    
    
    
    thresh_dirs_img=nib.Nifti2Image(thresh_dirs_array[0:counter,0:3,0],AFD_img.affine)
    
    
    nib.save(thresh_dirs_img, myargs.fixel_dir_name[0]+"/directions.nii") 
    
    
    t1_fixel_img=nib.Nifti2Image(t1_fixel_array[0:counter,0,0],AFD_img.affine) 

    
    nib.save(t1_fixel_img, myargs.fixel_dir_name[0]+"/t1_fixel.nii") 

    print('*********')
    print('All done.')
    print('*********')

        
#Visualize: 
#the visualization is done in voxel space == world space for no angulation and isotropic steps and therefore only looks correct for isotropic voxels (and a right handed coordinate system, which nifti is.)
#color code by T1
                
                
#DO: color code by AFD/(1/T1)? g-like?

#DO: remove underscores

if (myargs.sortT1): #sort T1s for other analysis:
    
    if (sagittal and xz): #x and z
        print("\n\nCrossing fiber T1s:\n")
        print("~x                   ~z\n")
        for j in range(len(voxels[0])):  
            if (number_of_fibers_array[voxels[0][j],voxels[1][j],voxels[2][j]]==2): #if two fibers
                if (abs(fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],0])>abs(fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],2])): #more along x than z       
                    ccT1=float(T1_array[voxels[0][j],voxels[1][j],voxels[2][j],0])
                    cingT1=float(T1_array[voxels[0][j],voxels[1][j],voxels[2][j],1])
                else:
                    ccT1=float(T1_array[voxels[0][j],voxels[1][j],voxels[2][j],1])
                    cingT1=float(T1_array[voxels[0][j],voxels[1][j],voxels[2][j],0])
                print("%f           %f" % (ccT1, cingT1))
        
        print("\n\nSingle fiber T1s:\n")
        print("~x                   \n")  
        for j in range(len(voxels[0])):         
            if (number_of_fibers_array[voxels[0][j],voxels[1][j],voxels[2][j]]==1): #single fiber      
                if (abs(fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],0])>abs(fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],2])): #x
                    print("%f" % float(T1_array[voxels[0][j],voxels[1][j],voxels[2][j],0]))
                
        
        print("\n~z                   \n")  
        for j in range(len(voxels[0])):         
            if (number_of_fibers_array[voxels[0][j],voxels[1][j],voxels[2][j]]==1): #single fiber      
                if (abs(fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],0])<=abs(fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],2])): #z
                    print("%f" % float(T1_array[voxels[0][j],voxels[1][j],voxels[2][j],0]))
                    
    if (axial and xz): #axial x&z             
        print("\n\nCrossing fiber T1s:\n")
        print("~x                   ~z\n")
        for j in range(len(voxels[0])):  
            if (number_of_fibers_array[voxels[0][j],voxels[1][j],voxels[2][j]]==2): #if two fibers
                if (abs(fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],0])>abs(fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],2])): #more along x than z      
                    ccT1=float(T1_array[voxels[0][j],voxels[1][j],voxels[2][j],0])
                    cingT1=float(T1_array[voxels[0][j],voxels[1][j],voxels[2][j],1])
                else:
                    ccT1=float(T1_array[voxels[0][j],voxels[1][j],voxels[2][j],1])
                    cingT1=float(T1_array[voxels[0][j],voxels[1][j],voxels[2][j],0])
                print("%f           %f" % (ccT1, cingT1))
        
        print("\n\nSingle fiber T1s:\n")
        print("~x                   \n")  
        for j in range(len(voxels[0])):         
            if (number_of_fibers_array[voxels[0][j],voxels[1][j],voxels[2][j]]==1): #single fiber      
                if (abs(fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],0])>abs(fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],2])): #x
                    print("%f" % float(T1_array[voxels[0][j],voxels[1][j],voxels[2][j],0]))
                
        
        print("\n~z                   \n")  
        for j in range(len(voxels[0])):         
            if (number_of_fibers_array[voxels[0][j],voxels[1][j],voxels[2][j]]==1): #single fiber      
                if (abs(fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],0])<=abs(fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],2])): #z
                    print("%f" % float(T1_array[voxels[0][j],voxels[1][j],voxels[2][j],0]))

    if (axial and xy): #axial x&y             
        print("\n\nCrossing fiber T1s:\n")
        print("~x                   ~y\n")
        for j in range(len(voxels[0])):  
            if (number_of_fibers_array[voxels[0][j],voxels[1][j],voxels[2][j]]==2): #if two fibers
                if (abs(fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],0])>abs(fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],1])): #more along x than y      
                    ccT1=float(T1_array[voxels[0][j],voxels[1][j],voxels[2][j],0])
                    cingT1=float(T1_array[voxels[0][j],voxels[1][j],voxels[2][j],1])
                else:
                    ccT1=float(T1_array[voxels[0][j],voxels[1][j],voxels[2][j],1])
                    cingT1=float(T1_array[voxels[0][j],voxels[1][j],voxels[2][j],0])
                print("%f           %f" % (ccT1, cingT1))
        
        print("\n\nSingle fiber T1s:\n")
        print("~x                   \n")  
        for j in range(len(voxels[0])):         
            if (number_of_fibers_array[voxels[0][j],voxels[1][j],voxels[2][j]]==1): #single fiber      
                if (abs(fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],0])>abs(fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],1])): #x
                    print("%f" % float(T1_array[voxels[0][j],voxels[1][j],voxels[2][j],0]))
                
        
        print("\n~y                   \n")  
        for j in range(len(voxels[0])):         
            if (number_of_fibers_array[voxels[0][j],voxels[1][j],voxels[2][j]]==1): #single fiber      
                if (abs(fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],0])<=abs(fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],1])): #y
                    print("%f" % float(T1_array[voxels[0][j],voxels[1][j],voxels[2][j],0]))



if (myargs.visualize):

    _renderwindow = vtk.vtkRenderWindow()
    _renderer = vtk.vtkRenderer()
    _renderwindowinteractor = vtk.vtkRenderWindowInteractor()
    _background_color = [1,1,1]        
    _renderer.SetBackground(_background_color[0],_background_color[1], _background_color[2])
    _renderwindow.AddRenderer(_renderer)
    _renderwindowinteractor.SetRenderWindow(_renderwindow)
    _renderwindowinteractor.Initialize()
    _interactor_style_trackball = vtk.vtkInteractorStyleTrackballCamera()
    _renderwindowinteractor.SetInteractorStyle(_interactor_style_trackball)
    _mapper = vtk.vtkPolyDataMapper()
    
    #now create the vectors to display:
    _polydata = vtk.vtkPolyData()
    _hedgehog = vtk.vtkHedgeHog()
    _points = vtk.vtkPoints()
    _scalars = vtk.vtkFloatArray()
    _vectors = vtk.vtkFloatArray()
    _lut = vtk.vtkLookupTable()
    #_lut.SetTableRange(vis_min, vis_max) #this isn't working
    _lut.Build()
    #print _lut.GetRange()
    _vectors.SetNumberOfComponents(3) 
    
    
        
    counter=0
    for j in range(len(voxels[0])):  
        

        
        #loop through number of fibers  
        for i in range(0,number_of_fibers_array[voxels[0][j],voxels[1][j],voxels[2][j]]):
            
            if sagittal:#for sagittal, 0 is y, 1 is z, and 2 is x:
                _points.InsertPoint(counter, (voxels[2][j],voxels[0][j],voxels[1][j])) 
            else:#for transverse:    
                _points.InsertPoint(counter, (voxels[0][j],voxels[1][j],voxels[2][j]))     
        
            thisfiber=0.4*fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],3*i:3*i+3]


            
            thist1=float(T1_array[voxels[0][j],voxels[1][j],voxels[2][j],i])       
            
            _vectors.InsertTuple(counter,tuple(thisfiber)) 
              
            
            #_scalars.InsertNextValue(thist1)
            _scalars.InsertNextValue((thist1-vis_min)/vis_range)
            
            counter+=1
            if sagittal:#for sagittal, 0 is y, 1 is z, and 2 is x:
                _points.InsertPoint(counter, (voxels[2][j],voxels[0][j],voxels[1][j])) 
            else:#for transverse:    
                _points.InsertPoint(counter, (voxels[0][j],voxels[1][j],voxels[2][j]))  
            
            #fiber pointing in other way, so centered:         
            revfiber=-1.0*thisfiber                   
            _vectors.InsertTuple(counter,tuple(revfiber))
            
            
            #_scalars.InsertNextValue(thist1)
            _scalars.InsertNextValue((thist1-vis_min)/vis_range)
            
            counter+=1
            
    
    #DO: implement window/level interaction
    
    _polydata.SetPoints(_points)
    _polydata.GetPointData().SetVectors(_vectors)
    _polydata.GetPointData().SetScalars(_scalars)
    
    _hedgehog.SetInput(_polydata)
    #hedgehog.SetScaleFactor(1)
    
    
    #if you want tubes instead of lines:
    #_tubefilter = vtk.vtkTubeFilter()
    #_tubefilter.SetInput(_hedgehog.GetOutput())  
    #_tubefilter.SetRadius(0.08)
    #_mapper.SetInput(_tubefilter.GetOutput())
    
    #if just lines (preferred):
    _mapper.SetInput(_hedgehog.GetOutput())
    
    
    
    _mapper.ScalarVisibilityOn()
    _mapper.SetLookupTable(_lut)
    _mapper.SetColorModeToMapScalars()
    _mapper.SetScalarModeToUsePointData()
    
    _actor = vtk.vtkActor()
    _actor.SetMapper(_mapper)
    _actor.GetProperty().SetLineWidth(10)
    _renderer.AddActor(_actor)
    
    #colourbar:
    _lut.SetHueRange(0.7,0.0)
    
    
    _colourbar = vtk.vtkScalarBarActor()
    _colourbar.SetLookupTable(_lut)
    _colourbar.SetNumberOfLabels(0)
    _renderer.AddActor2D(_colourbar)
    

    #show 3 orth planes of anatomical:
#this is totally unfinished, starting at HERE
#mrview seems to be OK for most of our needs
    
    #just read in with nibabel and put in vtk polydata as in mincdiffusion
#    underlay_img=nib.load("dcm2nii/20181009_130402MP2RAGE1mms003a1001_strides123.nii")
#    underlay_image_data_vtk=vtk.vtkImageData()
#    #HERE need mapping from corner of anatomical underlay to the corner of voxel space of T1
#    #same mapping for spacing (setting 0.5)
#    underlay_image_data_vtk.SetOrigin(0,0,0)
#    underlay_image_data_vtk.
#      _image_data_vtk->SetSpacing((double)image_data->GetXStep(),(double)image_data->GetYStep(),(double)image_data->GetZStep());
#      _image_data_vtk->SetDimensions(image_data->GetXSize(),image_data->GetYSize(),image_data->GetZSize());
#      _image_data_vtk->SetScalarTypeToFloat();//CLAUDE this fixes seg fault!!!!!!!
#    
#      // should keep a pointer of this for float_array->Delete();   // CLAUDE
#      vtkFloatArray *float_array=vtkFloatArray::New();
#      /*float_array->Allocate( image_data->GetXSize() * image_data->GetYSize() *
#		             image_data->GetZSize() );   // CLAUDE*/
#      
#
#      for(x=0;x<image_data->GetXSize();x++)
#	{
#	  for(y=0;y<image_data->GetYSize();y++) 
#	    { 
#	      for(z=0;z<image_data->GetZSize();z++) 
#		{
#	
#		  ijk[0]=x;
#		  ijk[1]=y;
#		  ijk[2]=z;
#
#		  float_array->InsertValue(_image_data_vtk->ComputePointId(ijk),image_data->GetValue(x,y,z));
#
#		  
#		}
#	    }
#	}
#      
#      _image_data_vtk->GetPointData()->SetScalars(float_array);
#    
# 
# 
#    picker = vtk.vtkCellPicker()
#   
#    
# 
#    image_plane_widgets_x=vtk.vtkImagePlaneWidget()
#    image_plane_widgets_x.SetInput(reader.GetOutput())
#    image_plane_widgets_x.DisplayTextOn()
#    image_plane_widgets_x.RestrictPlaneToVolumeOn()
#    image_plane_widgets_x.SetPlaneOrientationToXAxes()
#    image_plane_widgets_x.SetSliceIndex(reader.GetOutput().GetCenter()[0])
# 	  
#     #image_plane_widgets_x.SetLookupTable(add a lut)
# 	
#    image_plane_widgets_x.SetPicker(picker)
#     
#    image_plane_widgets_x.SetInteractor(_renderwindowinteractor)      
#    image_plane_widgets_x.On()
#    
#    
       
    _renderer.Render()

    _renderwindowinteractor.Start()
    

    
