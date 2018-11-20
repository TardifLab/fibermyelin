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

#constants:
EPSILON=9E-9

#hardcoded stuff:

#file orientation (if False, assumes axial)
global sagittal 
sagittal=False

global set_Dpar_equal
set_Dpar_equal=False #have to set in other file too

global bedpostx
bedpostx=False

global just_b0
just_b0=False #have to set in other file too

#for visualization W/L:
#asparagus 500-800,1500
#human brain 550-575-600,725-750-800
global vis_min
vis_min=650
global vis_max
vis_max=750
global vis_range
vis_range=vis_max-vis_min


#end hardcoded stuff



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

parser.add_argument('-vis', dest='visualize', help='visualize previously computed fiber T1s', required=False, action='store_true')
parser.add_argument('-sort', dest='sortT1', help='print out previously computed fiber T1s in 2-fiber voxels, sorted by orientation', required=False, action='store_true')
parser.add_argument('-t1', dest='t1_image_filename', help='T1 filename: output filename if doing full computation; pre-computed T1 filename if -vis option is selected', required=True, nargs=1)
parser.add_argument('-Dpar', dest='Dpar_image_filename', help='D_parallel filename', required=False, nargs=1)
parser.add_argument('-mask', dest='mask_image_filename', help='mask filename, for computation or visualization. If for visualization, must be the same mask previously used for computation', required=True, nargs=1)
parser.add_argument('-vic', dest='vic_image_filename', help='vic filename, for computation only.', required=False, nargs=1)
parser.add_argument('-fixel',dest='fixel_dir_name', help='T1 fixel directory name; currently implemented only for axial images!!', required=True, nargs=1)
#not doing this here, assume user took care of it, then we input fixel2voxel-ed version:
#parser.add_argument('-fixel', dest='fixel_dirname', help='fixel directory from mrtrix, must contain directions.mif, afd.mif, index.mif', required=True, nargs=1)              
parser.add_argument('-afd', dest='afd_image_filename', help='AFD filename, currently required even for -vis', required=True, nargs=1)
parser.add_argument('-afdthresh', dest='AFD_thresh', help='AFD threshold, required.  For -vis option, AFD threshold must be that used to compute the pre-computed T1 map', required=True, nargs=1)
parser.add_argument('-dirs', dest='dirs_image_filename', help='directions filename, currently required even for -vis', required=True, nargs=1)
parser.add_argument('-IRdiff', dest='IR_diff_image_filename', help='IR diffusion weighted image series; required for full computation',  required=False, nargs=1)                     
parser.add_argument('-TIs', dest='TIs_filename', help='text file of TIs for each slice; required for full computation',  required=False, nargs=1)                   
parser.add_argument('-bvals', dest='bvals_filename', help='text file of bvals; required for full computation',  required=False, nargs=1)                     
parser.add_argument('-bvecs', dest='bvecs_filename', help='text file of bvecs; required for full computation',  required=False, nargs=1)                     
parser.add_argument('-TR', dest='TR', help='nominal TR for short-TR steady state equation (ms)',required=False,nargs=1)


myargs = parser.parse_args()


#currently reads in .nii files
#assumes already done:
#standard registration? be careful because different slices have different TIs, etc.
#everything has to be sampled the same way, and like the IR-diff data, i.e., that number of slices. don't change the slicing because it determines TI
#i.e. check whether the AFD, etc. have the same smapling and fix them if not.

#whether or not we use the steady-state (non-infinite TR) equation is hardcoded (to True) in FiberT1Solver::IRDiffEqn\n\
#it requires TR to be input
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
    IR_diff_img=nib.load(myargs.IR_diff_image_filename[0])
    vic_img=nib.load(myargs.vic_image_filename[0])
    
#For visualize, sort, and full computation:  
AFD_img=nib.load(myargs.afd_image_filename[0])    
max_fibers=AFD_img.header.get_data_shape()[3]    



fiber_dirs_img=nib.load(myargs.dirs_image_filename[0]) 



#FSL BEDPOSTX version:

if (bedpostx):
    max_fibers=2
    AFD_img_fsl=[]
    fiber_dirs_img_fsl=[]
    AFD_img_fsl.append(nib.load("bedpost_outputs_strides/mean_f1samples_strides.nii"))
    AFD_img_fsl.append(nib.load("bedpost_outputs_strides/mean_f2samples_strides.nii"))
    fiber_dirs_img_fsl.append(nib.load("bedpost_outputs_strides/dyads1_strides.nii"))
    fiber_dirs_img_fsl.append(nib.load("bedpost_outputs_strides/dyads2_strides.nii"))
    
mask_img=nib.load(myargs.mask_image_filename[0])    
    
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
    
for j in range(len(voxels[0])):
       
    #get number of super-threshold AFDs=number_of_fibers            
    number_of_fibers=0        
    for l in range(max_fibers):                  
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
    
   
       
    for j in range(len(voxels[0])):

        number_of_fibers=number_of_fibers_array[voxels[0][j],voxels[1][j],voxels[2][j]]
        print('number_of_fibers %i:' % number_of_fibers)
        
        vic=vic_img.get_data()[voxels[0][j],voxels[1][j],voxels[2][j]]
        
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
            
        
      
        
        
        print(fiber_dirs)
        print("AFDs:")
        print(AFD)
        
        
        IR_DWIs=np.zeros((number_of_TIs, number_of_DWIs),float)
        
        for i in range(0,number_of_TIs):
            IR_DWIs=IR_diff_img.get_data()[voxels[0][j],voxels[1][j],voxels[2][j],:]
        
        #TIs is the TIs for this slice coordinate  **if axial, slices in z**          
        #**if sagittal, slices in x, which is still stored in dim 2**  
        
        
        TIs=TIs_allslices[voxels[2][j]]
                
                
        
        
        
        #TRs=TRs_allslices[voxels[2][j]] #HERE do this when have the whole TR history in signal comp.
 
       
        if (number_of_fibers>0):
            t1solver=FiberT1Solver()
            
            #to test fixing D, we have to giev it the voxel coord
            
            t1solver.SetInputData(fiber_dirs,AFD,IR_DWIs,TIs,grad_table,vic,TR,voxels[0][j],voxels[1][j],voxels[2][j],sagittal,set_Dpar_equal)
            T1s=np.zeros(number_of_fibers)
            Dparfit=np.zeros(number_of_fibers)
            #T1sandDparfit=np.zeros([number_of_fibers,2])
            T1sandDparfit=t1solver.GetT1s()
            
        #hack for just_b0:
        if (just_b0):
            number_of_fibers=1
            
        for i in range(0,number_of_fibers):
            if (T1sandDparfit!=None):
                if (set_Dpar_equal):
                    T1s[i]=T1sandDparfit[i]
                    Dparfit[i]=T1sandDparfit[number_of_fibers]
                    print('T1: %d' % T1s[i])
                    print('Dpar: %f' % Dparfit[i])
                else:
                    T1s[i]=T1sandDparfit[2*i]
                    Dparfit[i]=T1sandDparfit[2*i+1]
                    print('T1: %d' % T1s[i])
                    print('Dpar: %f' % Dparfit[i])
                T1_array[voxels[0][j],voxels[1][j],voxels[2][j],i]=T1s[i]
                
                Dparfit_array[voxels[0][j],voxels[1][j],voxels[2][j],i]=Dparfit[i]
        
        #optional: this plots the single voxel right now:
        #t1solver.VisualizeFibers()
    
    
    #output the T1 map.  can visualize like an afd file with mrtrix (DO: make sparse version, voxel2fixel)
    #start as a nonsparse copy of AFD file, DO: might want to create sparse file here
    T1_img = nib.Nifti1Image(T1_array, AFD_img.affine, AFD_img.header)
    nib.save(T1_img, myargs.t1_image_filename[0])

    Dparfit_img = nib.Nifti1Image(Dparfit_array, AFD_img.affine, AFD_img.header)
    if (myargs.Dpar_image_filename is None):
        nib.save(Dparfit_img, "Dpar.nii")
    else:    
        nib.save(Dparfit_img, myargs.Dpar_image_filename[0])
    
 
#maybe this until the output_fixel if should be inside that if, but it does no harm here   
        
T1_array_zeroed=np.zeros(AFD_img.header.get_data_shape()) 

for j in range(len(voxels[0])):  
    counter=0
    for l in range(max_fibers): 
                         
        if AFD_img.get_data()[voxels[0][j],voxels[1][j],voxels[2][j],l]>AFD_thresh:            
            T1_array_zeroed[voxels[0][j],voxels[1][j],voxels[2][j],l]=T1_array[voxels[0][j],voxels[1][j],voxels[2][j],counter]
            counter+=1
        else:
            T1_array_zeroed[voxels[0][j],voxels[1][j],voxels[2][j],l]=0.0
          
if (False):   
    #output the T1s to match the directions initially input:
    #the intent is to then use voxel2fixel, but it fails.  It writes the first T1 to all fixels.
    T1_img_zeroed = nib.Nifti1Image(T1_array_zeroed, AFD_img.affine, AFD_img.header)
    nib.save(T1_img_zeroed, "t1_zeroed.nii")   
     
    
    

   
#output the directions, t1, and index files in sparse format. 
#I'm doing this for the thresholded data
#index is the size of AFD, first frame
    
     
new_size=[AFD_img.header.get_data_shape()[0],AFD_img.header.get_data_shape()[1], AFD_img.header.get_data_shape()[2],2]

index_array=np.zeros(new_size)

#temporarily make these as enormous as possible:


thresh_dirs_array=np.zeros([AFD_img.header.get_data_shape()[0]*AFD_img.header.get_data_shape()[1]*AFD_img.header.get_data_shape()[2],3,1])

t1_fixel_array=np.zeros([AFD_img.header.get_data_shape()[0]*AFD_img.header.get_data_shape()[1]*AFD_img.header.get_data_shape()[2],1,1])

counter=0
#for i in range(AFD_img.header.get_data_shape()[0]*AFD_img.header.get_data_shape()[1]*AFD_img.header.get_data_shape()[2])
for i in range(len(voxels[0])):#for just the mask
    fibercounter=0
    for l in range(max_fibers):   
                    
        if AFD_img.get_data()[voxels[0][i],voxels[1][i],voxels[2][i],l]>AFD_thresh: 
            if (index_array[voxels[0][i],voxels[1][i],voxels[2][i],0]==0):#set on first fiber
                
                index_array[voxels[0][i],voxels[1][i],voxels[2][i],1]=counter
            thresh_dirs_array[counter,0,0]=fiber_dirs_img.get_data()[voxels[0][i],voxels[1][i],voxels[2][i],3*l]
            thresh_dirs_array[counter,1,0]=fiber_dirs_img.get_data()[voxels[0][i],voxels[1][i],voxels[2][i],3*l+1]
            thresh_dirs_array[counter,2,0]=fiber_dirs_img.get_data()[voxels[0][i],voxels[1][i],voxels[2][i],3*l+2]
            t1_fixel_array[counter,0,0]=T1_array_zeroed[voxels[0][i],voxels[1][i],voxels[2][i],l]
            
            counter+=1
            fibercounter+=1
            index_array[voxels[0][i],voxels[1][i],voxels[2][i],0]=fibercounter
            
     
if not os.path.exists(myargs.fixel_dir_name[0]):
    os.makedirs(myargs.fixel_dir_name[0])
           
index_img = nib.Nifti1Image(index_array,AFD_img.affine)


nib.save(index_img, myargs.fixel_dir_name[0]+"/index.nii")



thresh_dirs_img=nib.Nifti1Image(thresh_dirs_array[0:counter,0:3,0],AFD_img.affine)


nib.save(thresh_dirs_img, myargs.fixel_dir_name[0]+"/directions.nii") 


t1_fixel_img=nib.Nifti1Image(t1_fixel_array[0:counter,0,0],AFD_img.affine) 


nib.save(t1_fixel_img, myargs.fixel_dir_name[0]+"/t1_fixel.nii") 

    
#Visualize: 
#the visualization is done in voxel space == world space for no angulation and isotropic steps and therefore only looks correct for isotropic voxels (and a right handed coordinate system, which nifti is.)
#color code by T1
                
                
#DO: color code by AFD/(1/T1)? g-like?

#DO: remove underscores

if (myargs.sortT1): #sort T1s for other analysis:
    
    if (sagittal): #currently coded for x and z
        print("\n\nCrossing fiber T1s:\n")
        print("~x                   ~z\n")
        for j in range(len(voxels[0])):  
            if (number_of_fibers_array[voxels[0][j],voxels[1][j],voxels[2][j]]==2): #if two fibers
                if (fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],0]>fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],2]): #more along x than z       
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
                if (fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],0]>fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],2]): #x
                    print("%f" % float(T1_array[voxels[0][j],voxels[1][j],voxels[2][j],0]))
                
        
        print("\n~z                   \n")  
        for j in range(len(voxels[0])):         
            if (number_of_fibers_array[voxels[0][j],voxels[1][j],voxels[2][j]]==1): #single fiber      
                if (fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],0]<=fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],2]): #z
                    print("%f" % float(T1_array[voxels[0][j],voxels[1][j],voxels[2][j],0]))
                    
    else: #axial. currently coded for x&y             
        print("\n\nCrossing fiber T1s:\n")
        print("~x                   ~y\n")
        for j in range(len(voxels[0])):  
            if (number_of_fibers_array[voxels[0][j],voxels[1][j],voxels[2][j]]==2): #if two fibers
                if (fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],0]>fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],1]): #more along x than y      
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
                if (fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],0]>fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],1]): #x
                    print("%f" % float(T1_array[voxels[0][j],voxels[1][j],voxels[2][j],0]))
                
        
        print("\n~y                   \n")  
        for j in range(len(voxels[0])):         
            if (number_of_fibers_array[voxels[0][j],voxels[1][j],voxels[2][j]]==1): #single fiber      
                if (fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],0]<=fiber_dirs_array[voxels[0][j],voxels[1][j],voxels[2][j],1]): #z
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
    

    