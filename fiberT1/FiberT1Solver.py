# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 13:23:23 2017

@author: jcampbel
"""


import vtk
import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
from scipy import linalg
import random

#hardcoded global vars:
global plot_fit
plot_fit=False
global sagittal
sagittal=False #must also be set in calling script
global just_b0
just_b0=False
global fix_D_phantom3
fix_D_phantom3=False #if True, see below in code for hardcoded values for phantom3
#should get rid of the above and just change params D to args D when we fix D
global simulate
simulate=False
global fix_vic
fix_vic=True

#import scipy
#print(scipy.__version__)
from scipy.optimize import least_squares #least_squares #can't import? why? 

class FiberT1Solver:
    """class to compute T1 for each fiber"""
    
    def __init__(self):
        """this class doesn't require anything yet."""
        
    def GetT1s(self):
        
       
        
        
        #set initial estimates for the parameters:       
        self.init_params=np.zeros(2*self.number_of_fibers+2) #T1 and Dpar for each fiber        
        self.lowerbounds=np.zeros(2*self.number_of_fibers+2)
        self.upperbounds=np.zeros(2*self.number_of_fibers+2)
        for i in range(0,self.number_of_fibers):
            self.init_params[2*i]=700 #T1 in ms
            self.lowerbounds[2*i]=200#0
            self.upperbounds[2*i]=4000#np.inf#
            self.init_params[2*i+1]=1.5E-3 #Dpar in mm2/s            
            self.lowerbounds[2*i+1]=0.1E-3#0#
            self.upperbounds[2*i+1]=5.0E-3#np.inf
        
        
        #additional Johnson noise term: 
        
        self.init_params[2*self.number_of_fibers]=0.0001 #neta for Johnson noise term. could constrain to something reasonable based on data    
       
        #TEST! not using it:
        #self.lowerbounds[2*self.number_of_fibers]=0
        #self.upperbounds[2*self.number_of_fibers]=9E-9
        
        self.lowerbounds[2*self.number_of_fibers]=0 #we will square it, so reduce search
        self.upperbounds[2*self.number_of_fibers]=np.inf
        
        #unknown M0: this depends very much on the acquisition
        #use first signal point (first TI, b=0) and init T1, assume long TR
        self.init_params[2*self.number_of_fibers+1]=np.absolute(self.IR_DWIs[0]/(1-2*np.exp(-1.0*self.TIs[0]/700)))
        #self.init_params[2*self.number_of_fibers+1]=400 
        self.lowerbounds[2*self.number_of_fibers+1]=0
        self.upperbounds[2*self.number_of_fibers+1]=np.inf
        
        #we have number of TIs * number of bvalues observations, need to string them all out and add the constants to each observation
        #fastest varying will be diffusion, then TI
        args=np.zeros((self.number_of_TIs*self.number_of_diff_encodes,6+6*self.number_of_fibers))        
        #args[:,0]=bvals        
        #args[:,1]=observations
        #args[:,2]=TR
        #args[:,3]=not used right now
        #args[:,4]=TIs
        #args[:,5]=number of fibers
        #args[:,6]=vic
        #args[:,7,13,7+6*i,...]=AFD(fiber i - 0 offset)        
        #args[:,8,14,8+6*i,...]=gradient x direction in coordinate space with fiber dir along x
        #args[:,9,15,9+6*i,...]=gradient y direction in coordinate space with fiber dir along x
        #args[:,10,16,10+6*i,...]=gradient z direction in coordinate space with fiber dir along x
        #args[:,11,17,11+6*i,...]=Dpar(fiber) (used if fix_D_phantom3) 
        #args[:,12,18,12+6*i,...]=Dpar(fiber) (used if fix_D_phantom3)
        
        
        
        if (just_b0): #we are just going to repeat the b=0 images
            for i in range(0,self.number_of_TIs):
                #bvals are all zero so leave that as is (init to zero)
                #TIs:
                args[i*self.number_of_diff_encodes:(i+1)*self.number_of_diff_encodes,4]=np.ones(self.number_of_diff_encodes)*self.TIs[i] 
                
        else:
            for i in range(0,self.number_of_TIs):     
                #bvals:
                args[i*self.number_of_diff_encodes:(i+1)*self.number_of_diff_encodes,0]=self.grad_table.bvals[0:self.number_of_diff_encodes]
               
                #TIs:
                args[i*self.number_of_diff_encodes:(i+1)*self.number_of_diff_encodes,4]=np.ones(self.number_of_diff_encodes)*self.TIs[i] 
                #args[i*self.number_of_diff_encodes:(i+1)*self.number_of_diff_encodes,2]=np.ones(self.number_of_diff_encodes)*self.TRs[i] 
        
                
        #constants self.number_of_fibers,self.AFDs,self.fiber_dirs:
        args[:,5]=self.number_of_fibers*np.ones(self.number_of_TIs*self.number_of_diff_encodes)
        args[:,6]=self.vic*np.ones(self.number_of_TIs*self.number_of_diff_encodes)
        for i in range(0,self.number_of_fibers):
            args[:,7+6*i]=self.AFDs[i]*np.ones(self.number_of_TIs*self.number_of_diff_encodes)
            
        
        
        
        #string out the observations with diffusion fastest varying
        counter=0
        for i in range(0,self.number_of_diff_encodes):             
            for j in range(0,self.number_of_TIs):                   
                if (self.IR_DWIs[i]<9E-9):
                    print "IR DWI image value 0"                    
                    return None
                if (just_b0): # a bit of a hack, we will repeat:
                    args[j*self.number_of_diff_encodes+i,1]=self.IR_DWIs[j]
                    
                    
                else:                   
                    args[j*self.number_of_diff_encodes+i,1]=self.IR_DWIs[counter] 
                    
                counter+=1


        
        #just using the nominal TR:
        args[:,2]=self.TR*np.ones(self.number_of_TIs*self.number_of_diff_encodes)

        #create the g-vector for each fiber in coord system aligned along that fiber:        
        for k in range(0,self.number_of_fibers):  
            
            
            fx=self.fiber_dirs[k][0]            
            fy=self.fiber_dirs[k][1]            
            fz=self.fiber_dirs[k][2]
            
            
                            
            v_orth1=GetOrthVector([fx, fy,fz])
            v_orth2=np.cross([fx, fy, fz],v_orth1)
            
            
                           
            #transform g into coord system of (f,v_orth1,v_orth2):
                        
            xfm_matrix=linalg.inv(np.array([[fx, fy, fz], v_orth1, v_orth2]))
            
            if (fix_D_phantom3):
                #get Dperp and Dpar from tensor fit in single fiber voxels with appropriate orientation and AFD>0.6
                #for phantom3:            
                if (self.vox1>42):
                    if (fx>0.5):   #old fixed cord  
                        args[:,11+6*k]=0.769E-3*np.ones(self.number_of_TIs*self.number_of_diff_encodes) #Dpar
                        args[:,12+6*k]=0.403E-3*np.ones(self.number_of_TIs*self.number_of_diff_encodes) #Dperp
                    else: #cuprizone cord
                        args[:,11+6*k]=0.543E-3*np.ones(self.number_of_TIs*self.number_of_diff_encodes)
                        args[:,12+6*k]=0.252E-3*np.ones(self.number_of_TIs*self.number_of_diff_encodes)
                elif (self.vox1<=42 and self.vox1>34): 
                    if (fx>0.5): #new fixed cord:
                        args[:,11+6*k]=0.625E-3*np.ones(self.number_of_TIs*self.number_of_diff_encodes)
                        args[:,12+6*k]=0.256E-3*np.ones(self.number_of_TIs*self.number_of_diff_encodes) 
                    else: #cuprizone cord
                        args[:,11+6*k]=0.543E-3*np.ones(self.number_of_TIs*self.number_of_diff_encodes)
                        args[:,12+6*k]=0.252E-3*np.ones(self.number_of_TIs*self.number_of_diff_encodes)    
                else: 
                    if (fx>0.5): #fresh cord
                        args[:,11+6*k]=0.567E-3*np.ones(self.number_of_TIs*self.number_of_diff_encodes)
                        args[:,12+6*k]=0.215E-3*np.ones(self.number_of_TIs*self.number_of_diff_encodes)  
                    else: #cuprizone cord
                        args[:,11+6*k]=0.543E-3*np.ones(self.number_of_TIs*self.number_of_diff_encodes)
                        args[:,12+6*k]=0.252E-3*np.ones(self.number_of_TIs*self.number_of_diff_encodes)    
                
            
            #make gnew for each observation
            for j in range(0,self.number_of_diff_encodes):  
                                
                
                #bvecs are in **either voxel space or PRS**: same for phantom3 case  
                #we convert the axes here for non-transverse
                #this doesn't handle any angulation!!!
                #fiber directions are in world space, ==voxel space once we exchange the axes if no angulation
                if (sagittal):
                    gtest=self.grad_table.bvecs[j,0:3] 
                    g=np.ones(3)
                    g[0]=gtest[1]
                    g[1]=gtest[2]
                    g[2]=gtest[0]
                else:#axial
                    g=self.grad_table.bvecs[j,0:3] 
                
                gnew=np.zeros(3)            
                gnew[0]=g[0]*xfm_matrix[0,0]+g[1]*xfm_matrix[1,0]+g[2]*xfm_matrix[2,0]
                gnew[1]=g[0]*xfm_matrix[0,1]+g[1]*xfm_matrix[1,1]+g[2]*xfm_matrix[2,1]
                gnew[2]=g[0]*xfm_matrix[0,2]+g[1]*xfm_matrix[1,2]+g[2]*xfm_matrix[2,2]
                
                #check some things:
                #print("fiber %f %f %f v_orth1 %f %f %f v_orth2 %f %f %f\ng %f %f %f gnew %f %f %f det %f" % (fx, fy, fz, v_orth1[0], v_orth1[1], v_orth1[2], v_orth2[0], v_orth2[1], v_orth2[2], g[0], g[1], g[2], gnew[0], gnew[1], gnew[2], linalg.det(np.array([[fx, fy, fz],v_orth1, v_orth2]))))
                
                
                #DO: is g normalized?  yes.
                #gnew=gnew/linalg.norm(gnew)
            
                #string out:
                for i in range(0,self.number_of_TIs):
                    if (args[i*self.number_of_diff_encodes+j,0]>9E-9):                        
                        args[i*self.number_of_diff_encodes+j,8+6*k]=gnew[0]
                        args[i*self.number_of_diff_encodes+j,9+6*k]=gnew[1]
                        args[i*self.number_of_diff_encodes+j,10+6*k]=gnew[2]
                        
                    else: #set explicitly to zero again for b=0
                        args[i*self.number_of_diff_encodes+j,8+6*k]=0.0
                        args[i*self.number_of_diff_encodes+j,9+6*k]=0.0
                        args[i*self.number_of_diff_encodes+j,10+6*k]=0.0
            
        
        
        #TEST: repeat the first point (b=0, first TI) N more times:        
        if (True):
            N=np.int(0.1*self.number_of_TIs*self.number_of_diff_encodes)
            newargs=np.zeros((self.number_of_TIs*self.number_of_diff_encodes+N,6+6*self.number_of_fibers))
            for i in range(self.number_of_TIs*self.number_of_diff_encodes):
                newargs[i,:]=args[i,:]
            for i in range(self.number_of_TIs*self.number_of_diff_encodes,self.number_of_TIs*self.number_of_diff_encodes+N):
                newargs[i,:]=args[0,:]      
            
        #what about repeating all of the b=0s?
        if (False):
            numrep=10;
            newargs=np.zeros(((self.number_of_TIs)*(self.number_of_diff_encodes+numrep),6+6*self.number_of_fibers))    
            for i in range(self.number_of_TIs*self.number_of_diff_encodes):
                newargs[i,:]=args[i,:]
                
            counter=0
            for i in range(self.number_of_TIs*self.number_of_diff_encodes,self.number_of_TIs*(self.number_of_diff_encodes+numrep)):                       
                newargs[i,:]=args[(counter%self.number_of_TIs)*self.number_of_diff_encodes,:]  
                counter+=1
                
                
                
        if (simulate):#simulate data for these input fibers and AFDs: 
            #default params, except different T1s and no neta
            new_params=np.copy(self.init_params)
            new_params[self.number_of_fibers]=0.0
            new_params[0]=700             
            if (self.number_of_fibers>1):
                new_params[2]=700
            print new_params    
            sim_data=newargs[:,1]+IRDiffEqn(new_params,newargs)
            random.seed()
            #real only:
            #newargs[:,1]=np.absolute(sim_data+[random.gauss(0,25) for i in range(len(sim_data))]) 
            #two channels:
            newargs[:,1]=np.sqrt(np.square(sim_data+[random.gauss(0,25) for i in range(len(sim_data))])+np.square([random.gauss(0,25) for i in range(len(sim_data))]))
               
        #fit the equation: there are a lot of options here
        #bounds could be implemented for 'lm'
        
            
        #trust region reflective:        
        #res_lsq = least_squares(IRDiffEqn, self.init_params, method='trf', bounds=tuple([self.lowerbounds,self.upperbounds]),args=args)
        
        #trf, more b0s, 3-point jacobian computation:
        res_lsq = least_squares(IRDiffEqn, self.init_params, method='trf', bounds=tuple([self.lowerbounds,self.upperbounds]),args=newargs, jac='3-point')
        
        
        #lm: 
        #res_lsq = least_squares(IRDiffEqn, self.init_params, method='lm',  args=args)   
        
        #lm, more b0s:
        #res_lsq = least_squares(IRDiffEqn, self.init_params, method='lm',  args=newargs)
    
     
        
        print "fit status %i" %  res_lsq.status       
        if (not res_lsq.success):            
            return None
        
        #print all fitted parameters:
        print(res_lsq.x)
        
        
        if (plot_fit):#look at the data:             
            # Create a figure instance
            fig = plt.figure(1, figsize=(9, 6))
            # Create an axes instance
            ax = fig.add_subplot(111)
            # Create the plot
            
            #everything, diffusion fastest varying:
            #ax.scatter(range(self.number_of_TIs*self.number_of_diff_encodes),args[:,1] , s=2, alpha=0.4)
            
            #phantom3 and normal 30 dir acq
            #just the b=0, z y x grad orientations             
            if (True):
                
                
                plotdata=np.zeros([self.number_of_TIs,4])
                thisDWI=[0, 21, 5, 30]
                for i in range(self.number_of_TIs):
                    for j in range(4):
                        plotdata[i,j]=args[i*self.number_of_diff_encodes+thisDWI[j],1]
                ax.plot(range(self.number_of_TIs),plotdata[:,0], 'k--')
                ax.plot(range(self.number_of_TIs),plotdata[:,1], 'b--')
                ax.plot(range(self.number_of_TIs),plotdata[:,2], 'g--')
                ax.plot(range(self.number_of_TIs),plotdata[:,3], 'r--')
                ax.set_title('All TIs, b=0 (black), ~z (blue), ~y (green), ~x (red) gradient orientations', fontsize=18)
                
                #now set to predicted signal:
                
                pred_sig_res=IRDiffEqn(res_lsq.x,args)
                
                
                
                for i in range(self.number_of_TIs):                   
                    for j in range(4):
                        plotdata[i,j]=args[i*self.number_of_diff_encodes+thisDWI[j],1]+pred_sig_res[i*self.number_of_diff_encodes+thisDWI[j]]
                ax.plot(range(self.number_of_TIs),plotdata[:,0], 'k-')
                ax.plot(range(self.number_of_TIs),plotdata[:,1], 'b-')
                ax.plot(range(self.number_of_TIs),plotdata[:,2], 'g-')
                ax.plot(range(self.number_of_TIs),plotdata[:,3], 'r-')
                
                textstr='dashed: data\nsolid: fit'
                ax.text(0.05,0.95,textstr)
                
                #DO: predicted for a more reasonable Dpar: run all (fit) a second time


            #phantom2
            #1 is +x-z, 2 is -x-z, 3 is +y-z
            if (False):
                plotdata=np.zeros([self.number_of_TIs,4])
                thisDWI=[0, 1, 2, 3]
                #ax.scatter(range(self.number_of_TIs),self.IR_DWIs[thisDWI*self.number_of_TIs:(thisDWI+1)*self.number_of_TIs] , s=2, alpha=0.4)
                for i in range(self.number_of_TIs):
                    for j in range(4):
                        plotdata[i,j]=args[i*self.number_of_diff_encodes+thisDWI[j],1]
                ax.plot(range(self.number_of_TIs),plotdata[:,0], 'k--')
                ax.plot(range(self.number_of_TIs),plotdata[:,1], 'b--')
                ax.plot(range(self.number_of_TIs),plotdata[:,2], 'g--')
                ax.plot(range(self.number_of_TIs),plotdata[:,3], 'r--') 
                ax.set_title('All TIs, b=0 (black), +x-z (blue), -x-z (green), +y-z (red) gradient orientations', fontsize=18)
            
            
            #plt.xlabel('TI')
            #plt.ylabel('signal')
            
            plt.show()
        
#==============================================================================
#         for i in range(0,self.number_of_fibers):
#             if (not just_b0):
#                 self.T1s[i]=res_lsq.x[2*i]
#                 Dparfit[i]=res_lsq.x[2*i+1]
#             else:
#                 self.T1s[i]=res_lsq.x[0]
#==============================================================================
        
        #return the entire fit:
        return res_lsq.x
        
    #DO: potentially just give inputs upon initialization, or will we reuse? give inputs to GetT1s?    
    def SetInputData(self,fiber_dirs,AFDs,IR_DWIs,TIs,grad_table,vic,TR,vox0,vox1,vox2):
        
        
                
        self.AFDs = AFDs
        
        #each fiber dir has a vector, fiber_dirs is 2D
        self.fiber_dirs = fiber_dirs  
               
        #number_of_fibers=size of AFD array
        self.number_of_fibers=len(self.AFDs)
        
        self.IR_DWIs = IR_DWIs
        
        self.TIs = TIs
        
        self.number_of_TIs = len(self.TIs)
                        
        self.grad_table = grad_table
        
        self.number_of_diff_encodes = len(self.grad_table.bvals)/self.number_of_TIs
        
        self.vic=vic
        
        self.TR = TR
        
        self.vox0=vox0
        self.vox1=vox1
        self.vox2=vox2
        
        #DO: assert size of IR_DWIs checks out
        
        #world space to voxel space transform is unity if aquired with no angulation and axial
    
        #DO: check whether mrtrix outputs in voxel space or world space (could be either, but it is always xyz, even if sag acq)
        #its worldspace is the same as nii, which is not the same as dicom
    
    def VisualizeFibers(self):
        #make vtk objects of the fibers in the right orientation, cylinders
        #color code by T1
        #color code by AFD/(1/T1)? g-like (not linearly proportional.)
    
        #set up Display Window:
        #DO: remove underscores: make these either of scoep just this function (then nothing), or this class (then self.)
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
        _vectors.SetNumberOfComponents(3) #DO: size of T1s/AFDs/fibers
        
        #for a test, one point, at origin.  eventually put multiple points at world space coordinates
        #loop through number of fibers  
        for i in range(0,self.number_of_fibers):
            _points.InsertPoint(2*i, (0,0,0))             
            _vectors.InsertTuple(2*i,tuple(self.fiber_dirs[i]))
            _scalars.InsertNextValue(float(self.T1s[i])/1000)       
            #other way, so centered: DO: -1*vector()
            revfiber=[-1*self.fiber_dirs[i][0], -1*self.fiber_dirs[i][1], -1*self.fiber_dirs[i][2]]        
            _points.InsertPoint(2*i+1, (0,0,0))        
            _vectors.InsertTuple(2*i+1,tuple(revfiber))       
            _scalars.InsertNextValue(float(self.T1s[i])/1000)
        
        
        #DO: implement window/level interaction
        
        _polydata.SetPoints(_points)
        _polydata.GetPointData().SetVectors(_vectors)
        _polydata.GetPointData().SetScalars(_scalars)
        
        _hedgehog.SetInput(_polydata)
        #hedgehog.SetScaleFactor(1)
        
        #DO: scale to fit world space voxel? this assumes isotropic, and fits voxel space voxel
        
        
        _tubefilter = vtk.vtkTubeFilter()
        _tubefilter.SetInput(_hedgehog.GetOutput())  
        _tubefilter.SetRadius(0.04)
        
        _mapper.SetInput(_tubefilter.GetOutput())
        _mapper.ScalarVisibilityOn()
        _mapper.SetLookupTable(_lut)
        _mapper.SetColorModeToMapScalars()
        _mapper.SetScalarModeToUsePointData()
        
        _actor = vtk.vtkActor()
        _actor.SetMapper(_mapper)
        _renderer.AddActor(_actor)
        _renderer.Render()
        _renderwindowinteractor.Start()
        


def IRDiffEqn(params,*args): #equation for residuals; params is vector of the unknowns
 
    
    if (len(np.shape(args)) != 2):
        args=args[0][:][:]
 
    #DO: transfer params args into reasonably named variables for readability
    steady_state=True #note this doesn't work for our purposes, if TR is short (and is not necessary, if TR is long)      
    number_of_obs=np.shape(args)[0]
    eq=np.zeros(number_of_obs)
    sig=np.zeros(number_of_obs)
    out=np.zeros(number_of_obs)
    #normalize the AFD:
    sum_AFD=0       
    numfibers=int(args[0][5]) 
    
    vic=args[0][6] # currently gets hardcoded below
    norm_AFD=np.zeros(numfibers)
    
    for i in range(0,numfibers):
        sum_AFD+=args[0][7+6*i]
    for i in range(0,numfibers):    
        norm_AFD[i]=args[0][7+6*i]/sum_AFD
        
    
    for h in range(0,number_of_obs):  
        if (not just_b0):
            for i in range(0,numfibers):             
                term1=norm_AFD[i]
                
                if (steady_state):
                    term2=1-2*np.exp(-args[h][4]/params[2*i])+np.exp(-args[h][2]/params[2*i])
                else:
                    term2=1-2*np.exp(-args[h][4]/params[2*i])
                                      
                #TEST! 
                if (False):                     
                    if (i==0):
                        term2=1-2*np.exp(-args[h][4]/643)                      
                    if (i==1):
                        term2=1-2*np.exp(-args[h][4]/718)   
                
                #using vic assumes same tortuosity for all fibers. 
                #DO: should add a term that is a factor of T1 (params) to incorp myelin
                
                
                
                
                if (not fix_D_phantom3):
                    #dummy vic for now, DO: remove when really input vic from NODDI 
                    if (fix_vic):
                        vic=0.4
                    else:
                        vic=args[h][6]
                    #using vic assumes same tortuosity for all fibers. 
                    #averages the intra- (0) and extra-axonal tensors to get Dperp                 
                    Dperp_factor=(1-vic)**2 #squared because we also do weighted average with ic (Dperp==0) compartment                
                    Dpar=params[2*i+1]
                    
                    #TEST!!!                    
                    if (False):
                        #Dpar=1.5E-3 
                        if (i==0):
                            Dpar=0.0012
                        if (i==1):
                            Dpar=0.0012
                    
                    Dperp=Dperp_factor*Dpar 
                elif (fix_D_phantom3):                     
                    Dpar=args[h][11+6*i] 
                    Dperp=args[h][12+6*i] 
                      
                D=np.zeros([3,3])
                
               
                                               
                #D in coord system of fiber dir and orth vectors (f,v_orth1,v_orth2)
                                               
                                               
                D=np.array([[Dpar, 0, 0],[0, Dperp, 0],[0, 0, Dperp]])
                            
                
                gnew=[args[h][8+6*i], args[h][9+6*i], args[h][10+6*i]]
                #DO as numpy matrix dot product  (note not matrix mult)
                Dterm=0.0
                for j in range(0,3):
                    for k in range (0,3):
                        Dterm+=D[j][k]*gnew[j]*gnew[k]
                        
                #print("Dterm %f" % Dterm)                  
                          
                term3=np.exp(-args[h][0]*Dterm)
                #print("term1 %f term2 %f term3 %f" % (term1, term2, term3))
                
                eq[h]+=term1*term2*term3
        else: #just_b0
            i=0 #do just once for first fiber
            if (steady_state):
                term2=1-2*np.exp(-args[h][4]/params[2*i])+np.exp(-args[h][2]/params[2*i])
            else:
                term2=1-2*np.exp(-args[h][4]/params[2*i])
            eq[h]+=term2                          
        
        #take magnitude, mult by M0, and for low SNR images (currently b>300 or signal<100), add a Johnson noise term:        
        #DO: implement a smart way of determining which TIs we need this for        
        #if (args[h][0]>300):
        #if (args[h][1]<100):
        if (True): #noise term for ALL images
        #if (False): #no noise term
            sig[h]=params[2*numfibers+1]*sqrt(eq[h]**2+params[2*numfibers]**2)         
        else:#no noise term
            sig[h]=params[2*numfibers+1]*sqrt(eq[h]**2) 
        
        
        out[h]=sig[h]-args[h][1]  
        
    #if (numfibers>1):
        #print("Dpar: %f %f; T1: %f %f; SOS residuals: %f" % (params[1], params[3], params[0], params[2], np.sum(np.square(out))))
    
    
    return out


        
def iszero(x):
    epsilon=9E-9 #adjust if inappropriate
    return abs(x)<epsilon
        
def GetOrthVector(v):
    if iszero(v[0]) and iszero(v[1]):
        if  iszero(v[3]):
            raise ValueError('zero vector')
        else:
            return [0,1,0]
    else:
        len=np.sqrt(v[1]**2+v[0]**2)
        return [-v[1]/len,v[0]/len,0]
    

