# -*- coding: utf-8 -*-
"""
Created on Tursday Aug 1st 2019

@author: ileppe
"""



import vtk
import numpy as np
from math import sqrt, cos, pi


import matplotlib.pyplot as plt
from scipy import linalg
import random

#this code includes a bunch of tests, e.g. fixing Dpar, etc...
#some of these tests are controlled by these global variables, others are buried deeper in the code.

#hardcoded global vars:
global plot_fit
plot_fit=False
global just_b0
just_b0=False #have to set here and in calling script fiberMTmyelin_pipeline.py

# Simulations
# they have been set up 2 different ways, either varying MTR and AFDs or varying tensor shape.
# --sim1:varying MTR and AFDs
#       - need input MTR file (all different combos)
#       - need input AFD file (all different combos)
#       - read in the tensor shape from AD and RD
# --sim2: varying tensor shape
#       - need ADin (all different combos)
#       - need RDin (all different combos)
#       - fixed difference in MTR set here in code
#       - AFD fixed in input file (same for all combos)
#       - the actual assumed tensor has to be hard-coded
# --sim3: varying the tensor shape but also multiple observations with different AFDs
#       - fixed MTR
#       - need ADin (all different combos)
#       - need RDin (all different combos)
#       - need input AFD file (all different combos but 1 estimate for all observations)
#       - simulate with many different tensors, but set to 0.7 in fitting

global sim1
global sim2
global sim3

#only 1 can be true
sim1 = False
sim2 = False
sim3 = True

global simulate #currently only for Dpar eq, not just b0,

#don't use this without editing the simulation code: it is hardcoded to MTR=0.3 right now
#this will make fibers have MTR=0.3,0.4, ...
global sim_different_MTs

#simulation: will read in values instead of assuming they are all MTR=0.3 (need an additional input)
global sim_input_MTs

# note that we still need to specify AD and RD, they will be assumed in the fitting
global sim_input_tensor

# this will read each AFD combo but assume it's all the same measurement. i.e. like COMMIT, it will use all
# these different observations to make an estimate of the MTRs (despite the tensor being fixed and potentially wrong)
global sim_afd_combo
global AD_h
global RD_h

if sim1:
    simulate = True
    sim_different_MTs = False
    sim_input_MTs = True
    sim_input_tensor = False
    sim_afd_combo = False

if sim2:
    simulate = True
    sim_different_MTs = True
    sim_input_MTs = False
    sim_input_tensor = True
    sim_afd_combo = False
    AD_h = 14.9e-4
    RD_h = 3.8e-4
if sim3:
    simulate = True
    sim_different_MTs = True
    sim_input_MTs = False
    sim_input_tensor = True
    sim_afd_combo = True
    AD_h = 14.9e-4
    RD_h = 3.8e-4

global sim_noise_level
sim_noise_level=0 #S0 is hardcoded to 500 below, so 10 is 2%, 15 is 3%
#sim_noise_level=0
global sim_S0
sim_S0=500
global plot_res #this plots the residuals by iterating over solutions, only use this with 1 voxel!
plot_res=False #plot_fit needs to be True as well
global dr_afd #this plots the residuals by iterating over different 'doctored' afds, only use this with 1 voxel!
dr_afd=False #plot_fit needs to be True as well
global avg_tensor #use a fixed tensor for all fibers
avg_tensor = True
global linear_fit #remember to change in fiberMTmyelin_pipeline too
linear_fit = True



#if True, see hardcoded vic value in IRDiffEq function.
#0.4: lower than healthy, but a better fit near the noise floor
#I recommend actually computing it instead
#actual values are ~0.4-0.6 in human..
#0.4 seems a better fit in marm
global fix_vic 
fix_vic=False
global fixed_vic
fixed_vic=0.6

#this is the only option right now:
#global mean_field_tort
#mean_field_tort=True
#DO: make a variable "tort" that can have different cases

#this was a one-time experiment in rat spinal cord cuprizone phantom
#if True, see below in code for hardcoded values for phantom3
global fix_D_phantom3
fix_D_phantom3=False 

global number_of_contrasts
number_of_contrasts = 1 #for MTR

global set_Dpar_equal
set_Dpar_equal=True#HERE this has to be true right now. have to set here and in calling script

if (just_b0):
    set_Dpar_equal=False #has to be for code structure and Dpar is irrelevant




#import scipy
#print(scipy.__version__)
from scipy.optimize import least_squares #least_squares




class FiberMTSolver:
    """class to compute MT for each fiber"""

    def __init__(self):
        """this class doesn't require anything yet."""

    def GetMTs(self):

        #set initial estimates for the parameters:
        self.init_params = np.zeros(self.number_of_fibers + 3)  # MT of each fiber
        self.lowerbounds = np.zeros(self.number_of_fibers + 3)
        self.upperbounds = np.zeros(self.number_of_fibers + 3)

        for i in range(0,self.number_of_fibers):
            self.init_params[i]=self.mtrDW  # MTR (optionally input on command line, otherwise 0.3)
            self.lowerbounds[i]=0#0
            self.upperbounds[i]=0.8#np.inf#100%

        #this sets the fiber MTs to different values if desired
        if (sim2 or sim3):
            if (self.number_of_fibers>1):
                for i in range(0,self.number_of_fibers):
                    self.init_params[i]=0.3+i*0.1 #MTR
            #print "Nominal MTRs %f and %f" % (self.init_params[0], self.init_params[1])

        #this sets the MTRs based on an input file (for simulation)
        if (sim1):
            for i in range(0, self.number_of_fibers):
                self.init_params[i] = self.MTin[i]  # MTR from input file


        #Dpar in mm2/s  (all fibers the same)
        self.init_params[self.number_of_fibers]=1.7E-3 #Dpar in mm2/s
        self.lowerbounds[self.number_of_fibers]=0.1E-3#0#
        self.upperbounds[self.number_of_fibers]=5.0E-3#np.inf
        if avg_tensor: #the better way would be to re-write the code so it doesn't fit this param at all...
            self.init_params[self.number_of_fibers] = self.AD  # Dpar in mm2/s
            self.lowerbounds[self.number_of_fibers] = self.AD-1E-6# don't let it vary
            self.upperbounds[self.number_of_fibers] = self.AD+1E-6# don't let it vary

        #additional Johnson noise term neta:
        self.init_params[self.number_of_fibers+1]=0.0
        self.lowerbounds[self.number_of_fibers+1]=0
        self.upperbounds[self.number_of_fibers+1]=5

        #unknown S0: this depends very much on the acquisition
        #to initialize roughly, use MToff and b=0, assume long TR DO change to steady state eq
        self.init_params[self.number_of_fibers+2] = 250

        if (simulate):
            self.init_params[self.number_of_fibers+2]=sim_S0

        #S0:
        self.lowerbounds[self.number_of_fibers+2]=0
        self.upperbounds[self.number_of_fibers+2]=np.inf

        if (just_b0):
            self.init_params=np.zeros(3)
            self.lowerbounds=np.zeros(3)
            self.upperbounds=np.zeros(3)
            self.number_of_fibers=1
            self.init_params[0]=0.3#MTR
            self.lowerbounds[0]=0#0
            self.upperbounds[0]=1#np.inf# 100%

          #additional Johnson noise term neta:
            self.init_params[self.number_of_fibers]=0.0
            self.lowerbounds[self.number_of_fibers]=0
            self.upperbounds[self.number_of_fibers]=np.inf

            #unknown S0: this depends very much on the acquisition
            #to initialize roughly, use MToff and b=0, assume long TR DO change to steady state eq
            self.init_params[self.number_of_fibers+2]=self.MT_DWIs[0]

            #S0:
            self.lowerbounds[self.number_of_fibers+1]=0
            self.upperbounds[self.number_of_fibers+1]=np.inf

        
        #else:      #if Dpar for each fiber HERE this option is unfinished!!!; need to put fibers at the end so that other params keep their indices, for b0 option too
        #    self.init_params=np.zeros(2*self.number_of_fibers+2) #MTR and Dpar for each fiber        
        #    self.lowerbounds=np.zeros(2*self.number_of_fibers+2)
        #    self.upperbounds=np.zeros(2*self.number_of_fibers+2)
        #    for i in range(0,self.number_of_fibers):
        #        self.init_params[2*i]=0.3 #MTR
        #        self.lowerbounds[2*i]=0#0
        #        self.upperbounds[2*i]=1#np.inf#
        #        self.init_params[2*i+1]=1.7E-3 #Dpar in mm2/s
        #        self.lowerbounds[2*i+1]=0.1E-3#0#
        #        self.upperbounds[2*i+1]=5.0E-3#np.inf
            
  
         
         
            #additional Johnson noise term neta:         
            self.init_params[2*self.number_of_fibers]=0.0 #7.544243 is actual noise mean in asparagus phantom
         
            #neta:
            self.lowerbounds[2*self.number_of_fibers]=0 
            self.upperbounds[2*self.number_of_fibers]=np.inf
         
            #unknown S0: this depends very much on the acquisition
            #to initialize roughly, use MToff and b=0, assume long TR DO change to steady state eq
            self.init_params[self.number_of_fibers+2]=self.MT_DWIs[0]
    
            #S0:
            self.lowerbounds[2*self.number_of_fibers+1]=0
            self.upperbounds[2*self.number_of_fibers+1]=np.inf


        
        #we have number of contrasts * number of bvalues observations (for MTR we have 2 contrasts MTon and MToff)
        # need to string them all out and add the constants to each observation
        #fastest varying will be diffusion, then MT
        
        args=np.zeros((self.number_of_contrasts*self.number_of_diff_encodes,10+6*self.number_of_fibers))
        #args[:,0]=bvals        
        #args[:,1]=MToff,MTon
        #args[:,2]=AD
        #args[:,3]=RD
        #args[:,4]=MTws (text file of MT weights)
        #args[:,5]=number of fibers
        #args[:,6]=vic
        #args[:,7,13,7+6*i,...]=AFD(fiber i - 0 offset)        
        #args[:,8,14,8+6*i,...]=gradient x direction in coordinate space with fiber dir along x
        #args[:,9,15,9+6*i,...]=gradient y direction in coordinate space with fiber dir along x
        #args[:,10,16,10+6*i,...]=gradient z direction in coordinate space with fiber dir along x
        #args[:,11,17,11+6*i,...]=Dpar(fiber) (used if fix_D_phantom3) 
        #args[:,12,18,12+6*i,...]=Dpar(fiber) (used if fix_D_phantom3)
        #args[:,13+6*(self.number_of_fibers-1)]=viso
        #args[:,14+6*(self.number_of_fibers-1)]=B1
        #args[:,15+6*(self.number_of_fibers-1)]=not used
        
        if (just_b0): #we are just going to repeat the b=0 images
            for i in range(0,self.number_of_contrasts):
                #bvals are all zero so leave that as is (init to zero)
                args[i*self.number_of_diff_encodes:(i+1)*self.number_of_diff_encodes,4]=np.ones(self.number_of_diff_encodes)*1
                
        else:
            
            for i in range(0,self.number_of_contrasts):     
                #bvals:
                for j in range(0,self.number_of_diff_encodes):
                    #bvals:
                    args[i*self.number_of_diff_encodes+j,0]=self.grad_table.bvals[j]
                    #MTws:
                    args[i*self.number_of_diff_encodes+j,4]=self.MTws[j]
                    
    
        #constants:
        args[:,2]=self.AD*np.ones(self.number_of_contrasts*self.number_of_diff_encodes)
        args[:,3]=self.RD*np.ones(self.number_of_contrasts*self.number_of_diff_encodes)
        args[:,5]=self.number_of_fibers*np.ones(self.number_of_contrasts*self.number_of_diff_encodes)
        args[:,6]=self.vic*np.ones(self.number_of_contrasts*self.number_of_diff_encodes)
        args[:,13+6*(self.number_of_fibers-1)]=self.viso*np.ones(self.number_of_contrasts*self.number_of_diff_encodes)
        args[:,14+6*(self.number_of_fibers-1)]=self.B1*np.ones(self.number_of_contrasts*self.number_of_diff_encodes)
        for i in range(0,self.number_of_fibers):
            if (sim_afd_combo): #don't keep the AFD constant for each diff encode (sim3)
                # we want these 11 combos
                afd_combos = (np.linspace(0,1,11),np.linspace(1,0,11)) #[0,1],[0.1,0.9],etc
                # repeat each combo the number of directions (69 in this case) (number_of_contrasts is 1 for MTR
                # but we've strung together MTon and MToff so repeat twice)
                afd_input = np.append(np.repeat(afd_combos,69,axis=1),np.repeat(afd_combos,69,axis=1),axis=1)
                args[:, 7 + 6 * i] = afd_input[i,:]
            else: # the regular case where the AFD is fixed in each voxel
                args[:, 7 + 6 * i] = self.AFDs[i] * np.ones(self.number_of_contrasts * self.number_of_diff_encodes)


        
        
        #string out the observations with diffusion fastest varying, the data is already stored this way
        if (just_b0): #we will repeat (not necessary but allows the rest of this code to be used):
            for i in range(0,self.number_of_diff_encodes):             
                for j in range(0,number_of_contrasts):                   
                    if (self.MT_DWIs[i]<9E-9):
                        print "MT DWI image value 0"
                        return None
                
                    args[j*self.number_of_diff_encodes+i,1]=self.MT_DWIs[j*self.number_of_diff_encodes]
  
        else:     

            args[:,1]=self.MT_DWIs
                    
        if (linear_fit):
            ## linear fit
            # S = S0*fi*e(-bD)*(1-MTRi)
            # AX = B
            # A = S0*fi*e(-bD) (constant if we fix the tensor)
            # X = 1-MTRi
            # B = S
            ## S0 is the signal without diffusion and without MT, average from input
            S0_ind = np.nonzero(np.logical_and(self.grad_table.bvals==0,self.MTws==0))
            S0 = np.mean(self.MT_DWIs[S0_ind])

            ## Dterm constant for each fiber (for real data and sim1)
            if (sim_input_tensor): #hard-coded tensor, specify the actual tensor on the command line with AD_in and RD_in (this is for sim2 and sim3)
                Din = np.zeros([3, 3, 2]) #2 different tensors
                Din[:, :, 0] = np.array([[self.ADin[0], 0, 0], [0, self.RDin[0], 0], [0, 0, self.RDin[0]]])
                Din[:, :, 1] = np.array([[self.ADin[1], 0, 0], [0, self.RDin[1], 0], [0, 0, self.RDin[1]]])
                Ain = np.zeros((self.number_of_contrasts * self.number_of_diff_encodes, self.number_of_fibers))

            #same average tensor for both, specified by AD and RD on input
            D = np.zeros([3, 3])
            # D in coord system of fiber dir and orth vectors (f,v_orth1,v_orth2)
            D = np.array([[self.AD, 0, 0], [0, self.RD, 0], [0, 0, self.RD]])

            Dterm_array = np.zeros((self.number_of_contrasts * self.number_of_diff_encodes, self.number_of_fibers))

            A = np.zeros((self.number_of_contrasts * self.number_of_diff_encodes, self.number_of_fibers))
            X = np.zeros((self.number_of_fibers,1))
            B = self.MT_DWIs #the observations

        # normalize the AFD:
        if sim3:  # don't bother, already normalized
            norm_AFD = afd_input
        else:
            norm_AFD = np.zeros(self.number_of_fibers)
            sum_AFD = np.sum(AFD)
            for i in range(0, self.number_of_fibers):
                norm_AFD[i] = self.AFDs[i] / sum_AFD


        #the nominal TR and TE:
        #args[:,2]=self.TR*np.ones(self.number_of_contrasts*self.number_of_diff_encodes)
        #args[:,15+6*(self.number_of_fibers-1)]=self.TE*np.ones(self.number_of_contrasts*self.number_of_diff_encodes)

        #create the g-vector for each fiber in coord system aligned along that fiber:        
        for k in range(0,self.number_of_fibers):  

            f = np.zeros(3)
            
            f[0]=self.fiber_dirs[k,0]
            f[1]=self.fiber_dirs[k,1]
            f[2]=self.fiber_dirs[k,2]
                            
            v_orth1=GetOrthVector(f)
            v_orth2=np.cross(f,v_orth1)

            #transformation into coord system of (f,v_orth1,v_orth2):
            xfm_forward=np.transpose([f,v_orth1,v_orth2])
            xfm_matrix=linalg.inv(xfm_forward)          

            if (fix_D_phantom3):
                #get Dperp and Dpar from tensor fit in single fiber voxels with appropriate orientation and AFD>0.6
                #for phantom3:            
                if (self.vox1>42):
                    if (fx>0.5):   #old fixed cord  
                        args[:,11+6*k]=0.769E-3*np.ones(self.number_of_contrasts*self.number_of_diff_encodes) #Dpar
                        args[:,12+6*k]=0.403E-3*np.ones(self.number_of_contrasts*self.number_of_diff_encodes) #Dperp
                    else: #cuprizone cord
                        args[:,11+6*k]=0.543E-3*np.ones(self.number_of_contrasts*self.number_of_diff_encodes)
                        args[:,12+6*k]=0.252E-3*np.ones(self.number_of_contrasts*self.number_of_diff_encodes)
                elif (self.vox1<=42 and self.vox1>34): 
                    if (fx>0.5): #new fixed cord:
                        args[:,11+6*k]=0.625E-3*np.ones(self.number_of_contrasts*self.number_of_diff_encodes)
                        args[:,12+6*k]=0.256E-3*np.ones(self.number_of_contrasts*self.number_of_diff_encodes) 
                    else: #cuprizone cord
                        args[:,11+6*k]=0.543E-3*np.ones(self.number_of_contrasts*self.number_of_diff_encodes)
                        args[:,12+6*k]=0.252E-3*np.ones(self.number_of_contrasts*self.number_of_diff_encodes)    
                else: 
                    if (fx>0.5): #fresh cord
                        args[:,11+6*k]=0.567E-3*np.ones(self.number_of_contrasts*self.number_of_diff_encodes)
                        args[:,12+6*k]=0.215E-3*np.ones(self.number_of_contrasts*self.number_of_diff_encodes)  
                    else: #cuprizone cord
                        args[:,11+6*k]=0.543E-3*np.ones(self.number_of_contrasts*self.number_of_diff_encodes)
                        args[:,12+6*k]=0.252E-3*np.ones(self.number_of_contrasts*self.number_of_diff_encodes)    

            #make gnew for each observation
            for j in range(0,self.number_of_diff_encodes):  
                #bvecs are in **either voxel space or PRS**: same for phantom3 case
                #we convert the axes here for non-transverse
                #this doesn't handle any angulation!!!
                #fiber directions are in world space, ==voxel space once we exchange the axes and account for strides if no angulation
                if (self.sagittal):                        
                    gtest=self.grad_table.bvecs[j,0:3] 
                    g=np.ones(3)
                    
                    #Q why was this here? this would take us from xyz to sag yzx
                    #we think the fiber dirs are in xyz.
                    #g[0]=gtest[1]
                    #g[1]=gtest[2]
                    #g[2]=gtest[0]
                    
                    
                    #put it in xyz:
                    g[0]=gtest[2]
                    g[1]=gtest[0]
                    g[2]=gtest[1]
                    
                else:#axial
                    g=self.grad_table.bvecs[j,0:3] 
                    
                
                gnew=np.zeros(3)            
                
                
                gnew[0]=xfm_matrix[0,0]*g[0]+xfm_matrix[0,1]*g[1]+xfm_matrix[0,2]*g[2]
                gnew[1]=xfm_matrix[1,0]*g[0]+xfm_matrix[1,1]*g[1]+xfm_matrix[1,2]*g[2]
                gnew[2]=xfm_matrix[2,0]*g[0]+xfm_matrix[2,1]*g[1]+xfm_matrix[2,2]*g[2]
                
                #check some things:
                #print("fiber %f %f %f v_orth1 %f %f %f v_orth2 %f %f %f\ng %f %f %f gnew %f %f %f det %f" % (f[0], f[1], f[2], v_orth1[0], v_orth1[1], v_orth1[2], v_orth2[0], v_orth2[1], v_orth2[2], g[0], g[1], g[2], gnew[0], gnew[1], gnew[2], linalg.det(np.array([f,v_orth1, v_orth2]))))

                #g is normalized
                #gnew=gnew/linalg.norm(gnew)

                # string out:
                for i in range(0, self.number_of_contrasts):
                    if (args[i * self.number_of_diff_encodes + j, 0] > 9E-9):
                        args[i * self.number_of_diff_encodes + j, 8 + 6 * k] = gnew[0]
                        args[i * self.number_of_diff_encodes + j, 9 + 6 * k] = gnew[1]
                        args[i * self.number_of_diff_encodes + j, 10 + 6 * k] = gnew[2]

                    else:  # set explicitly to zero again for b=0
                        args[i * self.number_of_diff_encodes + j, 8 + 6 * k] = 0.0
                        args[i * self.number_of_diff_encodes + j, 9 + 6 * k] = 0.0
                        args[i * self.number_of_diff_encodes + j, 10 + 6 * k] = 0.0


                if (linear_fit):
                    # DO as numpy matrix dot product  (note not matrix mult)
                    Dterm = 0.0
                    DtermIn = np.zeros(2) # create Dterm for each fiber in sim2 (stick with 2 fibers)
                    for i in range(0, 3):
                        for h in range(0, 3):
                            Dterm += D[i, h] * gnew[i] * gnew[h]
                            if (sim2 or sim3):
                                DtermIn[k] += Din[i, h, k] * gnew[i] * gnew[h]


                    #print("Dterm %f %i %i" % (Dterm,k,j))
                    Dterm_array[j,k]=Dterm

                    # create A, the coefficients A = S0*fi*e(-bD)
                    if (simulate):
                        S0 = sim_S0

                    if (sim2):
                        Ain[j, k] = S0 * norm_AFD[k] * np.exp(-self.grad_table.bvals[j] * DtermIn[k])
                        A[j, k] = S0 * norm_AFD[k] * np.exp(-self.grad_table.bvals[j] * Dterm)
                    if (sim3):
                        Ain[j, k] = S0 * norm_AFD[k,j] * np.exp(-self.grad_table.bvals[j] * DtermIn[k])
                        A[j, k] = S0 * norm_AFD[k,j] * np.exp(-self.grad_table.bvals[j] * Dterm)
                    else: #regular case, not simulation
                        A[j, k] = S0 * norm_AFD[k] * np.exp(-self.grad_table.bvals[j] * Dterm)


        #to weight b=0 more (instead of actually acquiring more), repeat N more times:
        #don't add b=0s, there are already some in the acquisition
        newargs = args
              
                
        if (simulate):#simulate data for these input fibers and AFDs:
            #NOTE this is just for Dpar eq case
            #default params, except potentially different MTRs and no neta

            new_params=np.copy(self.init_params)
            new_params[self.number_of_fibers+1]=0.0 #neta

            sim_data=MTDiffEqn_signal(new_params,newargs)

            if (linear_fit):
                #S = S0*fi*e(-bD)*(1-MTRi)
                # the 1st half has mtw=0 and the 2nd half mtw=1
                for m in range(0, self.MTws.shape[0]):
                    if (sim2 or sim3): #if we've created A with 2 different tensors
                        sim_data[m] = np.matmul(Ain[m,], (1 - self.MTws[m] * new_params[0:self.number_of_fibers]))
                    else:
                        sim_data[m] = np.matmul(A[m,], (1-self.MTws[m]*new_params[0:self.number_of_fibers]) ) #this does the matrix multiplications

            #this will just be the magnitude difference:
            #sim_data=newargs[:,1]+IRDiffEqn(new_params,newargs) # this is the equation (sim) result
            random.seed()
            #real only: we don't want to do this
            #newargs[:,1]=np.absolute(sim_data+[random.gauss(0,sim_noise_level) for i in range(len(sim_data))])

            #add noise on two channels: we do want to do this
            newargs[:,1]=np.sqrt(np.square(sim_data+[random.gauss(0,sim_noise_level) for i in range(len(sim_data))])+np.square([random.gauss(0,sim_noise_level) for i in range(len(sim_data))]))

            # little test to see if I can make the data so noisy that the MTR is very variable
            if (0):
                fig = plt.figure(1, figsize=(9, 6))
                ax = fig.add_subplot(111)
                noise_level_test=[0,10,20,30]
                color = ['k-','r--','g--','b--']
                for s in range(0,len(noise_level_test)):
                    B=np.sqrt(np.square(sim_data+[random.gauss(0,noise_level_test[s]) for i in range(len(sim_data))])+np.square([random.gauss(0,noise_level_test[s]) for i in range(len(sim_data))]))
                    S_noMT = B[np.logical_and(self.grad_table.bvals != 0, self.MTws == 0)]
                    S_MT = B[np.logical_and(self.grad_table.bvals != 0, self.MTws == 1)]
                    mtr_dir = (S_noMT-S_MT)/S_noMT #MTR at every dir
                    A_MT = A[np.logical_and(self.grad_table.bvals != 0, self.MTws == 1)]  # only MT and no b=0s
                    B_MT = B[np.logical_and(self.grad_table.bvals != 0, self.MTws == 1)]  # only MT and no b=0s
                    A_MT = np.vstack([A_MT, np.tile(norm_AFD, (5, 1))])  # add 5times to increase the weight
                    B_MT = np.append(B_MT, (1 - np.mean(mtr_dir)) * np.ones(5))
                    X = np.linalg.lstsq(A_MT, B_MT)
                    ax.plot(mtr_dir, color[s], linewidth=2, label="noise level " + str(noise_level_test[s]) + " fit mtr "+str(1-X[0]))
                    ax.legend()
                #ax.legend()
                plt.ylim(-1, 1)
                plt.show()
            ##finish little test

            #write out the simulated data to a file
            if (linear_fit):
                B = newargs[:,1]
                # check whether the average diffusion-w MTR makes sense
                S_noMT = np.mean(B[np.logical_and(self.grad_table.bvals!=0, self.MTws==0)])#average signal over all directions; no MT; no b=0
                S_MT = np.mean(B[np.logical_and(self.grad_table.bvals != 0, self.MTws == 1)]) #average signal over all directions; with MT; no b=0
                self.mtrDW = (S_noMT-S_MT)/S_noMT
                #now is this same or similar to (1-MTR)=fi(1-MTRi) or MTR = 1-(fi(1-MTRi))
                mtrDW_sim = 1-(np.matmul(self.AFDs, (1 - new_params[0:self.number_of_fibers])))
                print "mtrDW = %f and mtrDW_sim = %f" % (self.mtrDW,mtrDW_sim)

        #fit the equation: there are a lot of options here; user can modify this call

        #bounds could be implemented for 'lm'
        
            
        #trust region reflective:        
        #res_lsq = least_squares(IRDiffEqn, self.init_params, method='trf', bounds=tuple([self.lowerbounds,self.upperbounds]),args=args)
        if (linear_fit):
            A_MT = A[np.logical_and(self.grad_table.bvals!=0, self.MTws==1)] #only MT and no b=0s
            B_MT = B[np.logical_and(self.grad_table.bvals!=0, self.MTws==1)] #only MT and no b=0s
            if (not sim3):
                #now add the extra constraint of the spherically averaged MTR [1-MTR = fi(1-MTRi)]
                A_MT = np.vstack([A_MT,np.tile(norm_AFD,(5,1))]) # add 5times to increase the weight
                B_MT = np.append(B_MT, (1-self.mtrDW)*np.ones(5))
            X= np.linalg.lstsq(A_MT,B_MT) #only use the MT weighted data so only the 2nd half
            result = (1-X[0]) #1-MTR
            if (any(result<0)): #if MTR<0, set to MTR of DW
                result = result.fill(self.mtrDW)

        else:
            #trf, more b0s, 3-point jacobian computation:
            res_lsq = least_squares(MTDiffEqn, self.init_params, method='trf', bounds=tuple([self.lowerbounds,self.upperbounds]),args=newargs, jac='3-point')

            #res_lsq = least_squares(MTDiffEqn, self.init_params, method='trf',bounds=tuple([self.lowerbounds, self.upperbounds]), args=newargs, jac='3-point', verbose=2)

        #lm:

            #res_lsq = least_squares(MTDiffEqn, self.init_params, method='lm',args=newargs)
        
            #lm, more b0s:
            #res_lsq = least_squares(IRDiffEqn, self.init_params, method='lm',  args=newargs)
    
     
        
            print "fit status %i" %  res_lsq.status
            if (not res_lsq.success or res_lsq.status==0 ):      #status=0 means it reached max iterations


                return None
        
            #print all fitted parameters:
            print(res_lsq.x)
                
        
            #print("%d SSE %f" % (self.TIs[0], res_lsq.cost))
            #print args[:,4] #TIs
            #print("vic %f" % self.vic)
            #print("viso %f" % self.viso)
        
        if (plot_fit):#look at the data:             
            # Create a figure instance
            fig = plt.figure(1, figsize=(9, 6))
            # Create an axes instance
            ax = fig.add_subplot(111)
            # Create the plot
            
            #everything, diffusion fastest varying:
            #ax.scatter(range(self.number_of_contrasts*self.number_of_diff_encodes),args[:,1] , s=2, alpha=0.4)
            
            #normal 30 dir acq
            #just the b=0, z y x grad orientations             
            if (True):
                
                
                plotdata=np.zeros([self.number_of_contrasts+1,4]) #IRL +1 because I strung the contrast after each other
                if (self.number_of_diff_encodes >= 65): #there might additional b=0s
                    thisDWI = [0, 59, 37, 35]#these are the index for b=0, and closest to be along z,y,z
                elif (self.number_of_diff_encodes >= 31):
                    thisDWI=[0, 21, 5, 30] #these are the index for b=0, and closest to be along z,y,z
                else: #(self.number_of_diff_encodes == 21):
                    thisDWI = [0, 18, 2, 1]
                for i in range(self.number_of_contrasts+1):
                    for j in range(4):
                        plotdata[i,j]=newargs[i*self.number_of_diff_encodes/2+thisDWI[j],1]
                ax.plot([0,1],plotdata[:,0], 'k--') #[0 1] is MTon,MToff
                ax.plot([0,1],plotdata[:,1], 'b--') #z
                ax.plot([0,1],plotdata[:,2], 'g--') #y
                ax.plot([0,1],plotdata[:,3], 'r--') #x
                
                ymax=500
                plt.ylim(0, ymax*1.1)
                #ax.set_title('All TIs, b=0 (black), ~z (blue), ~y (green), ~x (red) gradient orientations', fontsize=18)
                ax.set_title('Raw (--) and Fitted (-) Data',
                             fontsize=18)
                plt.xlabel('MT-weighting (ON-OFF)',fontsize=14)
                plt.ylabel('Signal Intensity',fontsize=14)
                
                #now set to predicted signal:
                if linear_fit:
                    pred_sig_res = np.matmul(A,result)
                else:
                    pred_sig_res=MTDiffEqn(res_lsq.x,newargs)
                
                
                
                for i in range(self.number_of_contrasts+1):
                    for j in range(4):
                        plotdata[i,j]=newargs[i*self.number_of_diff_encodes/2+thisDWI[j],1]+pred_sig_res[i*self.number_of_diff_encodes/2+thisDWI[j]]
                ax.plot([0,1],plotdata[:,0], 'k-', linewidth=2, label='b=0') #[0 1] is MToff,MTon
                ax.plot([0,1],plotdata[:,1], 'b-',linewidth=2, label='z-oriented')
                ax.plot([0,1],plotdata[:,2], 'g-',linewidth=2, label='y-oriented')
                ax.plot([0,1],plotdata[:,3], 'r-',linewidth=2, label='x-oriented')
                ax.legend()

                #these vary a bit with acquisition
                #textstr='dashed: data\nsolid: fit\n\nactual directions:\n~x=[-0.97958797216415, 0.17135678231716, 0.10509198158979]\n~y=[0.20307792723178, 0.94054549932479, -0.27227476239204]\n~z=[-0.20379328727722, 0.17156073451042, 0.96386468410492]'
                #ax.text(0.9,250,textstr)
                
                #DO: predicted for a more reasonable Dpar: run all (fit) a second time
                # write fitted MT to plot
                textstr =''
                for t in range(self.number_of_fibers):
                    if linear_fit:
                        textstr = textstr + (r'$MT_%d=%.2f (f_%d=%.1f)$' % (t + 1, result[t], t + 1, self.AFDs[t]))
                    else:
                        textstr = textstr+(r'$MT_%d=%.2f (f_%d=%.1f)$' % (t+1,res_lsq.x[t],t+1,self.AFDs[t]))
                    if (t!=self.number_of_fibers-1):
                        textstr =textstr + ('\n')
                props = dict(boxstyle='round', facecolor='grey', alpha=0.5)
                ax.text(0.3, 0.8, textstr, transform=ax.transAxes, fontsize=18,
                        verticalalignment='top', bbox=props)


            #phantom2
            #1 is +x-z, 2 is -x-z, 3 is +y-z
            if (False):
                plotdata=np.zeros([self.number_of_contrasts,4])
                thisDWI=[0, 1, 2, 3]
                #ax.scatter(range(self.number_of_contrasts),self.IR_DWIs[thisDWI*self.number_of_contrasts:(thisDWI+1)*self.number_of_contrasts] , s=2, alpha=0.4)
                for i in range(self.number_of_contrasts):
                    for j in range(4):
                        plotdata[i,j]=args[i*self.number_of_diff_encodes+thisDWI[j],1]
                ax.plot(range(self.number_of_contrasts),plotdata[:,0], 'k--')
                ax.plot(range(self.number_of_contrasts),plotdata[:,1], 'b--')
                ax.plot(range(self.number_of_contrasts),plotdata[:,2], 'g--')
                ax.plot(range(self.number_of_contrasts),plotdata[:,3], 'r--') 
                ax.set_title('All TIs, b=0 (black), +x-z (blue), -x-z (green), +y-z (red) gradient orientations', fontsize=18)
            
            
            #plt.xlabel('TI')
            #plt.ylabel('signal')
            
            plt.show()

            #plot solution space (this now only works for 2 fibers)
            if(plot_res):
                # Create a figure instance
                fig2 = plt.figure(2, figsize=(9, 6))
                # Create an axes instance
                ax = fig2.add_subplot(111)

                mt1=np.arange(0, 1, 0.05).tolist()
                mt2=np.arange(0, 1, 0.05).tolist()
                MSE = np.zeros([len(mt1), len(mt2)])
                for g in range(len(mt1)):
                    for h in range(len(mt2)):
                        res_lsq.x[0]=mt1[g]
                        res_lsq.x[1]=mt2[h]
                        pred_sig_res=MTDiffEqn(res_lsq.x,newargs)
                        summation=0

                        for p in range (len(pred_sig_res)):
                            squared_difference = pred_sig_res[p]**2
                            summation = summation + squared_difference
                        MSE[g,h] = summation/len(pred_sig_res)
                plt.imshow(MSE, interpolation='none')
                plt.xlabel('MT2', fontsize=14)
                plt.ylabel('MT1', fontsize=14)
                ax.set_title('MSE')
                plt.xticks(range(len(mt2)), mt2)
                plt.yticks(range(len(mt1)), mt1)
                plt.clim(100, 130)
                plt.colorbar(orientation='vertical')



            if(dr_afd): # what if we use ariticial AFDs, are the fits better? (only 2 fibers for now)
                # Create a figure instance
                fig2 = plt.figure(2, figsize=(9, 6))
                # Create an axes instance
                ax2 = fig2.add_subplot(111)

                afd1=np.arange(0.1, 1, 0.1).tolist()
                afd2=np.arange(0.1, 1, 0.1).tolist()
                MSE = np.zeros([len(afd1), len(afd2)])
                MT1 = np.zeros([len(afd1), len(afd2)]) #computed MT values
                MT2 = np.zeros([len(afd1), len(afd2)]) #computed MT values
                tempargs=newargs

                for g in range(len(afd1)):
                    for h in range(len(afd2)):
                        #args[:,7,13,7+6*i,...]=AFD(fiber i - 0 offset)   
                        tempargs[:,7]=afd1[g]
                        tempargs[:,13]=afd2[h]
                        #redo the fit
                        res_lsq_tmp = least_squares(MTDiffEqn, self.init_params, method='trf', bounds=tuple([self.lowerbounds,self.upperbounds]),args=tempargs, jac='3-point')
                        pred_sig_res=MTDiffEqn(res_lsq_tmp.x,tempargs)
                        summation=0
                        #compute the mean-squared-error
                        for p in range (len(pred_sig_res)):
                            squared_difference = pred_sig_res[p]**2
                            summation = summation + squared_difference
                        MSE[g,h] = summation/len(pred_sig_res)
                        MT1[g,h] = res_lsq_tmp.x[0]
                        MT2[g,h] = res_lsq_tmp.x[1]
                plt.imshow(MSE, interpolation='none')
                plt.xlabel('AFD2', fontsize=14)
                plt.ylabel('AFD1', fontsize=14)
                ax2.set_title('MSE')
                plt.xticks(range(len(afd2)), afd2)
                plt.yticks(range(len(afd1)), afd1)
                plt.clim(100, 150)
                plt.colorbar(orientation='vertical')

                ## now plot the MT values
                # Create a figure instance
                fig3 = plt.figure(3, figsize=(9, 6))
                # Create an axes instance
                ax3 = fig3.add_subplot(111)
                plt.imshow(MT1, interpolation='none')
                plt.xlabel('AFD2', fontsize=14)
                plt.ylabel('AFD1', fontsize=14)
                ax3.set_title('MT1')
                plt.xticks(range(len(afd2)), afd2)
                plt.yticks(range(len(afd1)), afd1)
                plt.colorbar(orientation='vertical')
                plt.clim(0, 0.5)

                # Create a figure instance
                fig4 = plt.figure(4, figsize=(9, 6))
                # Create an axes instance
                ax4 = fig4.add_subplot(111)
                plt.imshow(MT2, interpolation='none')
                plt.xlabel('AFD2', fontsize=14)
                plt.ylabel('AFD1', fontsize=14)
                ax4.set_title('MT2')
                plt.xticks(range(len(afd2)), afd2)
                plt.yticks(range(len(afd1)), afd1)
                plt.colorbar(orientation='vertical')
                plt.clim(0, 0.5)

                #print out the best fit
                print("Best fit\nMSE.min %f" % MSE.min())
                tmp=np.argwhere(MSE==MSE.min())
                print ("MT1 %f" % MT1[tmp[0,0], tmp[0,1]])
                print ("MT2 %f" % MT2[tmp[0,0], tmp[0,1]])
                print ("afd1 %f" % afd1[tmp[0,0]])
                print ("afd2 %f" % afd2[tmp[0,1]])



#         for i in range(0,self.number_of_fibers):
#             if (not just_b0):
#                 self.T1s[i]=res_lsq.x[2*i]
#                 Dparfit[i]=res_lsq.x[2*i+1]
#             else:
#                 self.T1s[i]=res_lsq.x[0]
#==============================================================================
        if linear_fit:
            return result
        else:
            # compute the MSE
            pred_sig_res=MTDiffEqn(res_lsq.x,newargs) #predicted signal residual
            summation = 0
            for p in range(len(pred_sig_res)):
                squared_difference = pred_sig_res[p] ** 2
                summation = summation + squared_difference
            MSE = summation / len(pred_sig_res)
            print("MSE %f" % MSE)
       
            #return the entire fit:
            return res_lsq
        
    #DO: potentially just give inputs upon initialization, or will we reuse? give inputs to GetT1s?    
    def SetInputData(self,fiber_dirs,AFDs,MT_DWIs,MTws,grad_table,AD,RD,mtrDW,vic,vox0,vox1,vox2,sag,Dpareq,viso,B1,MTin,ADin, RDin):
        self.AFDs = AFDs
        
        #each fiber dir has a vector, fiber_dirs is 2D
        self.fiber_dirs = fiber_dirs  
               
        #number_of_fibers=size of AFD array
        self.number_of_fibers=len(self.AFDs)
        
        self.MT_DWIs = MT_DWIs
        self.MTws = MTws

        self.number_of_contrasts = number_of_contrasts
                        
        self.grad_table = grad_table
        self.AD = AD
        self.RD = RD

        if (sim2 or sim3):
            self.ADin = ADin
            self.RDin = RDin
            print "sim AD %f %f" %(ADin[0], ADin[1])
            print "sim RD %f %f" %(RDin[0], RDin[1])

        self.mtrDW = mtrDW
        
        self.number_of_diff_encodes = len(self.grad_table.bvals)/self.number_of_contrasts
        
        self.vic=vic
        self.viso=viso
        self.B1=B1

        self.MTin=MTin #this is used for simulation when MT is used an an input sim_input_MTs=True

             
        self.vox0=vox0
        self.vox1=vox1
        self.vox2=vox2
        
        
        self.sagittal=sag
        
        #set_Dpar_equal=Dpareq #this doesn't work. it retains the value set outside of the methods here
        
        #DO: assert size of IR_DWIs checks out
        
        #world space to voxel space transform is unity if aquired with no angulation and axial
        #note the added issue of "strided" voxel space
        
        #DO: check whether mrtrix outputs in voxel space or world space (could be either, but it is always xyz, even if sag acq)
        #its worldspace is the same as nii, which is not the same as dicom
    



def MTDiffEqn(params,*args): #equation for residuals; params is vector of the unknowns
 
    #notes:
    #this always uses the steady state equation
    #using vic assumes same tortuosity for all fibers. 
    #DO: should add a term to tortuosity computation that is a factor of T1 (params) to incorp myelin
    # (although tortuousity model is wrong, period, so need DIAMOND or something else
              
    
    if (len(np.shape(args)) != 2):
        args=args[0][:][:]
    
 
    
    number_of_obs = np.shape(args)[0]

    #fixed params for CSF fraction:
    CSF_T1=2900 #2900 from Ilana's computation (agrees with Hutter);
    CSF_D=3.0E-3
    
    #put args in reasonably named variables
    #this takes longer but is more readable
    #I'm having trouble extracting the right thing from args, so I'm doing a hack:
    bvals = np.zeros(number_of_obs)
    obs  =np.zeros(number_of_obs)
    MTws=np.zeros(number_of_obs)
    
    temp_array=args[0]
    numfibers=int(temp_array[5]) #repeated constant
    if (fix_vic):
        vic=fixed_vic  #global variable
    else:
        vic=temp_array[6] #repeated constant
    viso=temp_array[13+6*(numfibers-1)] #repeated constant
    B1=temp_array[14+6*(numfibers-1)] #repeated constant

    if sim3: #as many different AFD combos as observations
        AFD = np.zeros([number_of_obs, numfibers])
    else:
        AFD = np.zeros(numfibers)
        for j in range(0, numfibers):
            AFD[j] = temp_array[7 + 6 * j]
            # AFD[j]=tmp[j]

    
    gnewx=np.zeros([number_of_obs,numfibers])
    gnewy=np.zeros([number_of_obs,numfibers])
    gnewz=np.zeros([number_of_obs,numfibers])
    
    
    for i in range(0,number_of_obs):
        temp_array=args[i]
        bvals[i]=temp_array[0]
        obs[i]=temp_array[1]
        MTws[i]=temp_array[4]

        for j in range(0,numfibers):
            gnewx[i,j]=temp_array[8+6*j]
            gnewy[i,j]=temp_array[9+6*j]
            gnewz[i,j]=temp_array[10+6*j]
            if sim3:
                AFD[i,j] = temp_array[7 + 6 * j]
  
    
    eq=np.zeros(number_of_obs)
    sig=np.zeros(number_of_obs)
    out=np.zeros(number_of_obs)
    # normalize the AFD:
    if sim3: #don't bother, already normalized
        norm_AFD = AFD
    else:
        norm_AFD = np.zeros(numfibers)
        sum_AFD = np.sum(AFD)
        for i in range(0, numfibers):
            norm_AFD[i] = AFD[i] / sum_AFD


    for h in range(0,number_of_obs):  
        if (not just_b0):
            for i in range(0,numfibers):
                if sim3:
                    term1 = norm_AFD[h, i]
                else:
                    term1 = norm_AFD[i]
                    # term1=AFD[i]


                #params[i] is MT-term for fiber i
                if (MTws[h]>=1): # if there is MT weighting
                    term2 = 1-params[i]
                else:
                    term2 =1 #no MT weigthing

                #else:  #not set_Dpar_equal not implemented right  now


                if (fix_D_phantom3):                     
                    Dpar=args[h,11+6*i] 
                    Dperp=args[h,12+6*i] 
                elif (avg_tensor):
                    Dpar = temp_array[2] #AD
                    Dperp = temp_array[3] #RD
                
                else: #(not fix_D_phantom3), i.e., everything else:
                   
                    
                    if(set_Dpar_equal):
                        Dpar=params[numfibers]
                    else: #Dpar for each fiber
                        if(i>0):
                            #params[numfibers+3+i], 1st is right after MT for each fiber, others are at end
                            Dpar=params[numfibers+2+i]
                        else:
                            Dpar = params[numfibers]


                    #this is the low-density mean-field tortuosity approximation, and is probably incorrect for realistically tight axonal packing
                    #hardcoded for b=1000; HERE change if b is not 1000!
                    Dperp=-np.log((1-vic)*np.exp(-(1-vic)*Dpar*1500)+vic)/1500
                    #Dperp = Dperp/2 #IRL reduce the anisotropy

                      
                D=np.zeros([3,3])
                
               
                                               
                #D in coord system of fiber dir and orth vectors (f,v_orth1,v_orth2)
                D=np.array([[Dpar, 0, 0],[0, Dperp, 0],[0, 0, Dperp]])

                #g in coordinate system of fiber i:
                
                gnew=[gnewx[h,i], gnewy[h,i], gnewz[h,i]]
                
                #DO as numpy matrix dot product  (note not matrix mult)
                Dterm=0.0
                for j in range(0,3):
                    for k in range (0,3):
                        Dterm+=D[j,k]*gnew[j]*gnew[k]
                        
                print("Dterm %f" % Dterm)
                          
                term3=np.exp(-bvals[h]*Dterm)

                #print("term1 %f term2 %f term3 %f" % (term1, term2, term3))
                
                eq[h]+=term1*term2*term3

            #now add CSF term:
            #eq[h]=(1-viso)*eq[h]+viso*(np.exp(-bvals[h]*CSF_D)*SteadyStateT1Recov(params[numfibers+3], B1, TIs[h], TR, TE, CSF_T1))


        else: #just_b0

            term2=SteadyStateT1Recov(params[numfibers+3], B1, TIs[h], TR, TE, params[0])

            #GE ver
            #term2=1-2*np.exp(-1.0*TIs[h]/params[0])+np.exp(-1.0*TR/params[0])

            eq[h]+=term2                          
        
        
        
        
        #mult by S0, and add Johnson noise term neta:
        #params[numfibers+2] is S0
        #params[numfibers+1] is noise term neta, currently added for ALL images

        #if (set_Dpar_equal):
        if True: #this should work for set_Dpar_equal or not
            #sig[h]=params[numfibers+2]*eq[h]+params[numfibers+1]
            sig[h] = params[numfibers+2]*eq[h] #IRL take out neta for now
        elif (just_b0):
            sig[h]=sqrt((params[numfibers+1]*eq[h])**2+params[numfibers]**2)
        #else:  #Dpar not equal, not implemented
        #    sig[h]=sqrt((params[2*numfibers+1]*eq[h])**2+params[2*numfibers]**2)
        out[h] = sig[h]-obs[h]
        
    #if (numfibers>1):
        #print("Dpar: %f; MT: %f %f; S0: %f SOS residuals: %f" % (params[2], params[0], params[1], params[4], np.sum(np.square(out))))
    #else:
        #print("Dpar: %f; MT: %f; S0: %f SOS residuals: %f" % (params[1], params[0], params[3], np.sum(np.square(out))))

    
    return out


def MTDiffEqn_signal(params, *args):   #equation for predicted signal; params is vector of the unknowns

    # notes:
    # this always uses the steady state equation
    # using vic assumes same tortuosity for all fibers.
    # DO: should add a term to tortuosity computation that is a factor of T1 (params) to incorp myelin
    # (although tortuousity model is wrong, period, so need DIAMOND or something else

    if (len(np.shape(args)) != 2):
        args = args[0][:][:]

    number_of_obs = np.shape(args)[0]

    # fixed params for CSF fraction:
    CSF_T1 = 2900  # 2900 from Ilana's computation (agrees with Hutter);
    CSF_D = 3.0E-3

    # put args in reasonably named variables
    # this takes longer but is more readable
    # I'm having trouble extracting the right thing from args, so I'm doing a hack:
    bvals = np.zeros(number_of_obs)
    obs = np.zeros(number_of_obs)
    MTws = np.zeros(number_of_obs)

    temp_array = args[0]
    if sim3: #hard-code to 2 fibers
        numfibers = 2
    else:
        numfibers = int(temp_array[5])  # repeated constant
    if (fix_vic):
        vic = fixed_vic  # global variable
    else:
        vic = temp_array[6]  # repeated constant
    viso = temp_array[13 + 6 * (numfibers - 1)]  # repeated constant
    B1 = temp_array[14 + 6 * (numfibers - 1)]  # repeated constant

    if sim3: #as many different AFD combos as observations
        AFD = np.zeros([number_of_obs, numfibers])
    else:
        AFD = np.zeros(numfibers)
        for j in range(0, numfibers):
            AFD[j] = temp_array[7 + 6 * j]
            # AFD[j]=tmp[j]

    gnewx = np.zeros([number_of_obs, numfibers])
    gnewy = np.zeros([number_of_obs, numfibers])
    gnewz = np.zeros([number_of_obs, numfibers])

    for i in range(0, number_of_obs):
        temp_array = args[i]
        bvals[i] = temp_array[0]
        obs[i] = temp_array[1]
        MTws[i] = temp_array[4]

        for j in range(0, numfibers):
            gnewx[i, j] = temp_array[8 + 6 * j]
            gnewy[i, j] = temp_array[9 + 6 * j]
            gnewz[i, j] = temp_array[10 + 6 * j]
            if sim3:
                AFD[i,j] = temp_array[7 + 6 * j]

    eq = np.zeros(number_of_obs)
    sig = np.zeros(number_of_obs)
    out = np.zeros(number_of_obs)

    # normalize the AFD:
    if sim3: #don't bother, already normalized
        norm_AFD = AFD
    else:
        norm_AFD = np.zeros(numfibers)
        sum_AFD = np.sum(AFD)
        for i in range(0, numfibers):
            norm_AFD[i] = AFD[i] / sum_AFD


    for h in range(0, number_of_obs):
        if (not just_b0):
            for i in range(0, numfibers):
                if sim3:
                    term1 = norm_AFD[h,i]
                else:
                    term1 = norm_AFD[i]
                    # term1=AFD[i]

                if (set_Dpar_equal):

                    # params[i] is MT-term for fiber i
                    if (MTws[h] >= 1):  # if there is MT weighting
                        term2 = 1 - params[i]
                    else:
                        term2 = 1  # no MT weigthing

                # else:  #not set_Dpar_equal not implemented right  now

                if (fix_D_phantom3):
                    Dpar = args[h, 11 + 6 * i]
                    Dperp = args[h, 12 + 6 * i]
                elif (avg_tensor):
                    Dpar = args[h,2] #AD
                    Dperp = args[h,3] #RD

                else:  # (not fix_D_phantom3), i.e., everything else:

                    if (set_Dpar_equal):
                        Dpar = params[numfibers]
                    else:  # not implemented right now
                        # params[numfibers+3 + i] is Dpar2, Dpar3 etc
                        Dpar = params[numfibers+3 + i]

                    # this is the low-density mean-field tortuosity approximation, and is probably incorrect for realistically tight axonal packing
                    # hardcoded for b=1000; HERE change if b is not 1000!
                    Dperp = -np.log((1 - vic) * np.exp(-(1 - vic) * Dpar * 1500) + vic) / 1500
                    # Dperp = Dperp/2 #IRL reduce the anisotropy

                D = np.zeros([3, 3])

                # D in coord system of fiber dir and orth vectors (f,v_orth1,v_orth2)
                D = np.array([[Dpar, 0, 0], [0, Dperp, 0], [0, 0, Dperp]])

                # g in coordinate system of fiber i:

                gnew = [gnewx[h, i], gnewy[h, i], gnewz[h, i]]

                # DO as numpy matrix dot product  (note not matrix mult)
                Dterm = 0.0
                for j in range(0, 3):
                    for k in range(0, 3):
                        Dterm += D[j, k] * gnew[j] * gnew[k]

                # print("Dterm %f" % Dterm)

                term3 = np.exp(-bvals[h] * Dterm)

                # print("term1 %f term2 %f term3 %f" % (term1, term2, term3))

                eq[h] += term1 * term2 * term3

            # now add CSF term:
            # eq[h]=(1-viso)*eq[h]+viso*(np.exp(-bvals[h]*CSF_D)*SteadyStateT1Recov(params[numfibers+3], B1, TIs[h], TR, TE, CSF_T1))


        else:  # just_b0

            term2 = SteadyStateT1Recov(params[numfibers + 3], B1, TIs[h], TR, TE, params[0])

            # GE ver
            # term2=1-2*np.exp(-1.0*TIs[h]/params[0])+np.exp(-1.0*TR/params[0])

            eq[h] += term2

            # mult by S0, and add Johnson noise term neta:
        # params[numfibers+2] is S0
        # params[numfibers+1] is noise term neta, currently added for ALL images

        if (set_Dpar_equal):
            # sig[h]=params[numfibers+2]*eq[h]+params[numfibers+1]
            sig[h] = params[numfibers + 2] * eq[h]  # IRL take out neta for now
        elif (just_b0):
            sig[h] = sqrt((params[numfibers + 1] * eq[h]) ** 2 + params[numfibers] ** 2)
        else:  # Dpar not equal, not implemented
            sig[h] = sqrt((params[2 * numfibers + 1] * eq[h]) ** 2 + params[2 * numfibers] ** 2)
        out[h] = sig[h]

    # if (numfibers>1):
    # print("Dpar: %f %f; T1: %f %f; SOS residuals: %f" % (params[1], params[3], params[0], params[2], np.sum(np.square(out))))


    return out

def SignedIRDiffEqn(params,*args): #equation for predicted signal; params is vector of the unknowns
    #copied from IRDiffEqn(params,*args) with the final section modified to output signed signal not magnitude residual




    #notes:
    #this always uses the steady state equation
    #using vic assumes same tortuosity for all fibers.
    #DO: should add a term to tortuosity computation that is a factor of T1 (params) to incorp myelin
    # (although tortuousity model is wrong, period, so need DIAMOND or something else


    if (len(np.shape(args)) != 2):
        args=args[0][:][:]



    number_of_obs=np.shape(args)[0]

    #fixed params for CSF fraction:
    CSF_T1=2900 #2900 from Ilana's computation (agrees with Hutter);
    CSF_D=3.0E-3



    #put args in reasonably named variables
    #this takes longer but is more readable
    #I'm having trouble extracting the right thing from args, so I'm doing a hack:
    bvals=np.zeros(number_of_obs)
    obs=np.zeros(number_of_obs)
    TIs=np.zeros(number_of_obs)

    temp_array=args[0]
    TR=temp_array[2]#currently a constant; we don't have a sequence with a different but constant TR per slice. Need Bloch sim if it varies.
    numfibers=int(temp_array[5]) #repeated constant
    if (fix_vic):
        vic=fixed_vic  #global variable
    else:
        vic=temp_array[6] #repeated constant
    viso=temp_array[13+6*(numfibers-1)] #repeated constant
    B1=temp_array[14+6*(numfibers-1)] #repeated constant
    TE=temp_array[15+6*(numfibers-1)]

    AFD=np.zeros(numfibers)
    for j in range(0,numfibers):
        AFD[j]=temp_array[7+6*j]

    gnewx=np.zeros([number_of_obs,numfibers])
    gnewy=np.zeros([number_of_obs,numfibers])
    gnewz=np.zeros([number_of_obs,numfibers])


    for i in range(0,number_of_obs):
        temp_array=args[i]
        bvals[i]=temp_array[0]
        obs[i]=temp_array[1]
        TIs[i]=temp_array[4]
        for j in range(0,numfibers):
            gnewx[i,j]=temp_array[8+6*j]
            gnewy[i,j]=temp_array[9+6*j]
            gnewz[i,j]=temp_array[10+6*j]




    eq=np.zeros(number_of_obs)
    sig=np.zeros(number_of_obs)
    out=np.zeros(number_of_obs)

    #normalize the AFD:

    norm_AFD=np.zeros(numfibers)
    sum_AFD=np.sum(AFD)
    for i in range(0,numfibers):
        norm_AFD[i]=AFD[i]/sum_AFD


    for h in range(0,number_of_obs):
        if (not just_b0):
            for i in range(0,numfibers):
                term1=norm_AFD[i]

                if (set_Dpar_equal):
                    #params[i] is T1 for fiber i
                    #inversion efficiency IE is params[numfibers+3]

                    term2=SteadyStateT1Recov(params[numfibers+3], B1, TIs[h], TR, TE, params[i])
                    

                    #GE ver:
                    #term2=1-2*params[numfibers+3]*np.exp(-1.0*TIs[h]/params[i])+np.exp(-1.0*TR/params[i])



                #else:  #not set_Dpar_equal not implemented right  now
                #params[2*i] is T1 for fiber i



                    #GE ver:
                    #term2=1-2*np.exp(-1.0*TIs[h]/params[2*i])+np.exp(-1.0*TR/params[2*i])



                if (fix_D_phantom3):
                    Dpar=args[h,11+6*i]
                    Dperp=args[h,12+6*i]


                else: #(not fix_D_phantom3), i.e., everything else:


                    if (set_Dpar_equal):
                        Dpar=params[numfibers]
                    else:  #not implemented right now
                    #params[2*i+1] is Dpar
                        Dpar=params[2*i+1]





                    #this is the low-density mean-field tortuosity approximation, and is probably incorrect for realistically tight axonal packing
                    #hardcoded for b=1000; HERE change if b is not 1000!
                    Dperp=-np.log((1-vic)*np.exp(-(1-vic)*Dpar*1000)+vic)/1000



                D=np.zeros([3,3])



                #D in coord system of fiber dir and orth vectors (f,v_orth1,v_orth2)
                D=np.array([[Dpar, 0, 0],[0, Dperp, 0],[0, 0, Dperp]])




                #g in coordinate system of fiber i:

                gnew=[gnewx[h,i], gnewy[h,i], gnewz[h,i]]

                #DO as numpy matrix dot product  (note not matrix mult)
                Dterm=0.0
                for j in range(0,3):
                    for k in range (0,3):
                        Dterm+=D[j,k]*gnew[j]*gnew[k]

                #print("Dterm %f" % Dterm)

                term3=np.exp(-bvals[h]*Dterm)



                #print("term1 %f term2 %f term3 %f" % (term1, term2, term3))

                eq[h]+=term1*term2*term3

            #now add CSF term:
            eq[h]=(1-viso)*eq[h]+viso*(np.exp(-bvals[h]*CSF_D)*SteadyStateT1Recov(params[numfibers+3], B1, TIs[h], TR, TE, CSF_T1))


        else: #just_b0

            term2=SteadyStateT1Recov(params[numfibers+3], B1, TIs[h], TR, TE, params[0])

            #GE ver
            #term2=1-2*np.exp(-1.0*TIs[h]/params[0])+np.exp(-1.0*TR/params[0])

            eq[h]+=term2




        #take magnitude, mult by S0, and add Johnson noise term neta:
        #params[2*numfibers+1] is S0
        #params[2*numfibers] is noise term neta, currently added for ALL images


        if (set_Dpar_equal):
            sig[h]=params[numfibers+2]*eq[h]
        elif (just_b0):
            sig[h]=params[numfibers+1]*eq[h]
        #else:  #Dpar not equal, not implemented

        out[h]=sig[h]

    #if (numfibers>1):
        #print("Dpar: %f %f; T1: %f %f; SOS residuals: %f" % (params[1], params[3], params[0], params[2], np.sum(np.square(out))))


    return out

def SteadyStateT1Recov(IE, B1, TI, TR, TE, T1):
    #RF pulses
    theta1=IE*pi
    theta2=B1*pi/2
    theta3=B1*pi


    term1=(1-cos(theta1)*cos(theta3)*np.exp(-TR/T1)-cos(theta1)*(1-cos(theta3))*np.exp(-(TR-TE/2)/T1))/(1-cos(theta1)*cos(theta2)*cos(theta3)*np.exp(-TR/T1))

    term2=-1.0*(1-cos(theta1))/(1-cos(theta1)*cos(theta2)*cos(theta3)*np.exp(-TR/T1))

    term3=term1+term2*np.exp(-TI/T1)


    return term3
        
def iszero(x):
    epsilon=9E-9 #adjust if inappropriate
    return abs(x)<epsilon
        
def GetOrthVector(v):
    if iszero(v[0]) and iszero(v[1]):
        if  iszero(v[2]):# IRL I don't why v[3] and not v[2]?
            raise ValueError('zero vector')
        else:
            return [0,1,0]    
    else:
        len=np.sqrt(v[1]**2+v[0]**2)
        return [-1.0*v[1]/len,v[0]/len,0]
    


def VisualizeDirs(dirs): #dirs is 1D, 3*number of dirs, xyz fastest varying
    

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
    _vectors.SetNumberOfComponents(3) 
    
    #one point, at origin.  
    #loop through number of dirs  
    for i in range(0,len(dirs)/3):
        _points.InsertPoint(i, (0,0,0))             
        _vectors.InsertTuple(i,tuple([dirs[3*i], dirs[1+3*i], dirs[2+3*i]]))
        _scalars.InsertNextValue(i)
    
    
    _polydata.SetPoints(_points)
    _polydata.GetPointData().SetVectors(_vectors)
    _polydata.GetPointData().SetScalars(_scalars)
    
    _hedgehog.SetInput(_polydata)
    #hedgehog.SetScaleFactor(1)
    

    
    _mapper.SetInput(_hedgehog.GetOutput())
    _mapper.ScalarVisibilityOn()
    _mapper.SetLookupTable(_lut)
    _mapper.SetColorModeToMapScalars()
    _mapper.SetScalarModeToUsePointData()
    
    _actor = vtk.vtkActor()
    _actor.SetMapper(_mapper)
    _actor.GetProperty().SetLineWidth(6)
    _renderer.AddActor(_actor)
    _renderer.Render()
    _renderwindowinteractor.Start()
    