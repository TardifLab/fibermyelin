#!/usr/bin/python
#/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 13:43:18 2017

@author: jcampbel
"""

'''This is a pipeline to process IR-diffusion data'''

#test

#from pyminc.volumes.factory import *
from fiberT1.FiberT1Solver import FiberT1Solver
#import os
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
sagittal=True

#for visualization W/L:
global vis_min
vis_min=0
global vis_range
vis_range=2000

#long TR: TRa not used:
TR=8000

#end hardcoded stuff

print('This a script to process IR-diffusion data\n\
NB: it does not currently handle data acquired with angulation\n\
NB: it assumes the following has already been done:\n\
    -registration of all images and sampling to same voxel space\n\
    -AFD computation with MRtrix\n\
    -optionally, tensor and/or NODDI fit\n\
NB: vic is currently hardcoded to 0.4 (lower than healthy, but a better fit near the noise floor).\n\
    You must input vic file regardless:\n\
    just input your mask. Edit FiberT1Solver global var fix_vic to actually use vic.\n\
NB: parallel diffusivity is currently output in file Dpar.nii\n\
NB: there is currently a Johnson noise term added to all images: we need to think about this implementation.\n\
NB: the data are assumed to have been acquired with diffusion fastest varying, and unshuffled to TI fastest varying.\n\
    i.e., bvals/bvecs have diffusion fastest varying, and input IRdiff images do not.\n\
NB: the visualization W/L is a WIP\n\
    todo: "voxel2fixel" the T1 output to view it in mrview.\n\
NB: initial parameter values are hardcoded in FiberT1Solver::GetT1s,\n\
    as is the fitting algorithm: currently trust region reflective (\'trf\').\n\
    Upper and lower bounds are also hardcoded\n\
    Comment/uncomment appropriate lines to change fitting algorithm and fitting options.\n\
NB: the first TI, b=0 image is currently weighted significantly more than the rest.\n\
NB: future additions that aren\'t here: CSF compartment, GM compartment, short TR\n\
    (Bloch sim, using TRa file for spin history by slice)')


parser = argparse.ArgumentParser(description='')

parser.add_argument('-vis', dest='visualize', help='visualize previously computed fiber T1s', required=False, action='store_true')
parser.add_argument('-t1', dest='t1_image_filename', help='T1 filename: output filename if doing full computation; pre-computed T1 filename if -vis option is selected', required=True, nargs=1)
#parser.add_argument('-fixel', dest='fixel_dirname', help='fixel directory from mrtrix, must contain directions.mif, afd.mif, index.mif', required=True, nargs=1)                   
parser.add_argument('-mask', dest='mask_image_filename', help='mask filename, for computation or visualization. If for visualization, must be the same mask previously used for computation', required=True, nargs=1)
parser.add_argument('-vic', dest='vic_image_filename', help='vic filename, for computation only.', required=False, nargs=1)
#not doing this here, assume user took care of it:
#parser.add_argument('-fixel', dest='fixel_dirname', help='fixel directory from mrtrix, must contain directions.mif, afd.mif, index.mif', required=True, nargs=1)                   
parser.add_argument('-afd', dest='afd_image_filename', help='AFD filename, currently required even for -vis', required=True, nargs=1)
parser.add_argument('-afdthresh', dest='AFD_thresh', help='AFD threshold, required.  For -vis option, AFD threshold must be that used to compute the pre-computed T1 map', required=True, nargs=1)
parser.add_argument('-dirs', dest='dirs_image_filename', help='directions filename, currently required even for -vis', required=True, nargs=1)
parser.add_argument('-IRdiff', dest='IR_diff_image_filename', help='IR diffusion weighted image series; required for full computation',  required=False, nargs=1)                     
parser.add_argument('-TIs', dest='TIs_filename', help='text file of TIs for each slice; required for full computation',  required=False, nargs=1)                   
parser.add_argument('-bvals', dest='bvals_filename', help='text file of bvals; required for full computation',  required=False, nargs=1)                     
parser.add_argument('-bvecs', dest='bvecs_filename', help='text file of bvecs; required for full computation',  required=False, nargs=1)                     



myargs = parser.parse_args()


#currently reads in .nii files
#assumes already done:
#standard registration, probably with FSL eddy of *all* images, including the diffusion only dataset
#everything has to be sampled the same way, and like the IR-diff data, i.e., that number of slices. don't change the slicing because it determines TI

#whether or not we use the steady-state (non-infinite TR) equation is hardcoded (to True) in FiberT1Solver::IRDiffEqn\n\
#however, it won't work if slices shuffled and TR short enough to matter, bc Tr will vary too much for a slice: need Bloch sim


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

#do same to mask as to files that match it (first 3 dims)

#same for vic 
#for axial (the only place I have tested vic), I needed to set strides of the NODDI output to -1,2,3 





if myargs.visualize:
    T1_array=nib.load(myargs.t1_image_filename[0]).get_data()
else:
    IR_diff_img=nib.load(myargs.IR_diff_image_filename[0])
    
#For both visualize and full computation:  
AFD_img=nib.load(myargs.afd_image_filename[0])    
max_fibers=AFD_img.header.get_data_shape()[3]    

fiber_dirs_img=nib.load(myargs.dirs_image_filename[0]) 
    
mask_img=nib.load(myargs.mask_image_filename[0])    
    
voxels=np.where(mask_img.get_data()>=EPSILON)  

number_of_fibers_array=np.zeros(mask_img.header.get_data_shape(),int)
    
#make a new AFD and fiber_dirs array for thresholded fiber dirs (they aren't sequential,at least beyond first)
AFD_array=np.zeros(AFD_img.header.get_data_shape())
fiber_dirs_array=np.zeros(fiber_dirs_img.header.get_data_shape())

AFD_thresh=float(myargs.AFD_thresh[0])

vic_img=nib.load(myargs.vic_image_filename[0])

    
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

if not myargs.visualize:    
    
  
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
            
            t1solver.SetInputData(fiber_dirs,AFD,IR_DWIs,TIs,grad_table,vic,TR,voxels[0][j],voxels[1][j],voxels[2][j],sagittal)
            T1s=np.zeros(number_of_fibers)
            Dparfit=np.zeros(number_of_fibers)
            #T1sandDparfit=np.zeros([number_of_fibers,2])
            T1sandDparfit=t1solver.GetT1s()
            
        
        for i in range(0,number_of_fibers):
            if (T1sandDparfit!=None):
                T1s[i]=T1sandDparfit[2*i]
                Dparfit[i]=T1sandDparfit[2*i+1]
                print('T1: %d' % T1s[i])
                print('Dpar: %f' % Dparfit[i])
                T1_array[voxels[0][j],voxels[1][j],voxels[2][j],i]=T1s[i]
                Dparfit_array[voxels[0][j],voxels[1][j],voxels[2][j],i]=Dparfit[i]
        
        #optional: this plots the single voxel right now:
        #t1solver.VisualizeFibers()
    
    
    #output the T1 map.  can visualize like an afd file with mrtrix
    #start as a nonsparse copy of AFD file, DO: might want to create sparse file here
    T1_img = nib.Nifti1Image(T1_array, AFD_img.affine, AFD_img.header)
    nib.save(T1_img, myargs.t1_image_filename[0])

    Dparfit_img = nib.Nifti1Image(Dparfit_array, AFD_img.affine, AFD_img.header)
    nib.save(Dparfit_img, "Dpar.nii")
    
#Visualize: currently always done
#the visualization is done in voxel space == world space for no angulation and isotropic steps and therefore only looks correct for isotropic voxels (and a right handed coordinate system, which nifti is.)
#color code by T1
                
                
#DO: color code by AFD/(1/T1)? g-like?

#DO: remove underscores



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
_lut.SetRange(0.0, 1.0)
_lut.Build()
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
          
        _scalars.InsertNextValue((float(T1_array[voxels[0][j],voxels[1][j],voxels[2][j],i])-vis_min)/vis_range)
        #_scalars.InsertNextValue(thist1/1000)
        #to vis afd: don't scale:
        #_scalars.InsertNextValue(thist1)        
        counter+=1
        if sagittal:#for sagittal, 0 is y, 1 is z, and 2 is x:
            _points.InsertPoint(counter, (voxels[2][j],voxels[0][j],voxels[1][j])) 
        else:#for transverse:    
            _points.InsertPoint(counter, (voxels[0][j],voxels[1][j],voxels[2][j]))  
        
        #fiber pointing in other way, so centered:         
        revfiber=-1.0*thisfiber                   
        _vectors.InsertTuple(counter,tuple(revfiber))
        #HERE manually playing with W/L     
        _scalars.InsertNextValue((float(T1_array[voxels[0][j],voxels[1][j],voxels[2][j],i])-vis_min)/vis_range)
        #_scalars.InsertNextValue(thist1/1000)
        #to vis afd: don't scale:
        #_scalars.InsertNextValue(thist1)
        counter+=1
        

#DO: implement window/level interaction

_polydata.SetPoints(_points)
_polydata.GetPointData().SetVectors(_vectors)
_polydata.GetPointData().SetScalars(_scalars)

_hedgehog.SetInput(_polydata)
#hedgehog.SetScaleFactor(1)


#if you want tubes instead of lines:
_tubefilter = vtk.vtkTubeFilter()
_tubefilter.SetInput(_hedgehog.GetOutput())  
_tubefilter.SetRadius(0.05)
_mapper.SetInput(_tubefilter.GetOutput())

#if lines:
#_mapper.SetInput(_hedgehog.GetOutput())

_mapper.ScalarVisibilityOn()
_mapper.SetLookupTable(_lut)
_mapper.SetColorModeToMapScalars()
_mapper.SetScalarModeToUsePointData()

_actor = vtk.vtkActor()
_actor.SetMapper(_mapper)
_renderer.AddActor(_actor)
_renderer.Render()
_renderwindowinteractor.Start()
