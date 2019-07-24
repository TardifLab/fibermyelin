# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 13:23:23 2017

@author: jcampbel
"""



import vtk
import numpy as np
from math import sqrt, cos, pi


import matplotlib.pyplot as plt
from scipy import linalg
import random

#this code includes a bunch of tests, e.g. fixing Dpar, etc...
#some of these tests are controlled by these global variables, others are buried deeper in the code.

#hardcoded T1 and D for CSF compartment are in IRDiffEqn

#hardcoded global vars:

global plot_fit
plot_fit=False

global just_b0
just_b0=True #have to set here and in calling script fibermyelin_pipeline.py, a

global simulate #currently only for Dpar eq, not just b0,
simulate=False
global sim_noise_level
sim_noise_level=15 #S0 is hardcoded to 500 below, so 10 is 2%, 15 is 3%
global sim_S0
sim_S0=500

global set_noIE
set_noIE=True

#don't use this without editing the simulation code: it is hardcoded to T1=750 right now
#this will make fibers have T1=700,800, ...

global sim_different_T1s
sim_different_T1s=False

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



global set_Dpar_equal
set_Dpar_equal=True #HERE this has to be true right now because the code hasn't been kept up to date for al cases. have to set here and in calling script



#import scipy
#print(scipy.__version__)
from scipy.optimize import least_squares #least_squares


if (just_b0):
    set_Dpar_equal=False #for if strucuture to work; DO: clean this up

class FiberT1Solver:
    """class to compute T1 for each fiber"""

    def __init__(self):
        """this class doesn't require anything yet."""

    def GetT1s(self):

       #set initial estimates for the parameters:



        if (set_Dpar_equal):

            if (not set_noIE):
                self.init_params=np.zeros(self.number_of_fibers+4) #T1 for each fiber
                self.lowerbounds=np.zeros(self.number_of_fibers+4)
                self.upperbounds=np.zeros(self.number_of_fibers+4)
            else:
                self.init_params = np.zeros(self.number_of_fibers + 3)  # T1 for each fiber
                self.lowerbounds = np.zeros(self.number_of_fibers + 3)
                self.upperbounds = np.zeros(self.number_of_fibers + 3)

            for i in range(0,self.number_of_fibers):
               self.init_params[i]=750 #T1 in ms
               self.lowerbounds[i]=200#0
               self.upperbounds[i]=4000#np.inf#

           #this sets the fiber T1s to different values if desired
            if (simulate and sim_different_T1s):
               if (self.number_of_fibers>1):
                    for i in range(0,self.number_of_fibers):
                        self.init_params[i]=700+i*100 #T1 in ms


            #Dpar in mm2/s  (all fibers the same)
            self.init_params[self.number_of_fibers]=1.7E-3 #Dpar in mm2/s
            self.lowerbounds[self.number_of_fibers]=0.1E-3#0#
            self.upperbounds[self.number_of_fibers]=5.0E-3#np.inf

            #additional Johnson noise term neta:
            self.init_params[self.number_of_fibers+1]=0.0
            self.lowerbounds[self.number_of_fibers+1]=0
            self.upperbounds[self.number_of_fibers+1]=np.inf

            #unknown S0: this depends very much on the acquisition
            #to initialize roughly, use first signal point (first TI, b=0) and init T1, assume long TR DO change to steady state eq
            self.init_params[self.number_of_fibers+2]=np.absolute(self.IR_DWIs[0]/(1-2*np.exp(-1.0*self.TIs[0]/750)))


            if (simulate):
               self.init_params[self.number_of_fibers+2]=sim_S0

            #S0:
            self.lowerbounds[self.number_of_fibers+2]=0
            self.upperbounds[self.number_of_fibers+2]=np.inf

            if (not set_noIE):
            #Inversion Efficiency IE
                self.init_params[self.number_of_fibers+3]=1.0
                self.lowerbounds[self.number_of_fibers+3]=0.5 #assuming at least 90!
                self.upperbounds[self.number_of_fibers+3]=1.0

        elif (just_b0):
            self.init_params=np.zeros(3) #1 T1
            self.lowerbounds=np.zeros(3)
            self.upperbounds=np.zeros(3)

            self.number_of_fibers=1
            self.init_params[0]=750#T1 in ms
            self.lowerbounds[0]=200#0
            self.upperbounds[0]=4000#np.inf#

            #additional Johnson noise term neta:
            self.init_params[1]=0.0
            self.lowerbounds[1]=0
            self.upperbounds[1]=np.inf

            #unknown S0: this depends very much on the acquisition
            #use first signal point (first TI, b=0) and init T1, assume long TR DO change to steady state eq
            self.init_params[2]=np.absolute(self.IR_DWIs[0]/(1-2*np.exp(-1.0*self.TIs[0]/750)))

            if (simulate):
               self.init_params[2]=sim_S0

            #S0:
            self.lowerbounds[2]=0
            self.upperbounds[2]=np.inf

            if (not set_noIE):
                #Inversion Efficiency IE
                self.init_params[3]=1.0
                self.lowerbounds[3]=0.5 #assuming at least 90!
                self.upperbounds[3]=1.0

        else:      #if Dpar for each fiber HERE this option is unfinished!!!; need to put fibers at the end so that other params keep their indices, for b0 option too
            self.init_params=np.zeros(2*self.number_of_fibers+2) #T1 and Dpar for each fiber
            self.lowerbounds=np.zeros(2*self.number_of_fibers+2)
            self.upperbounds=np.zeros(2*self.number_of_fibers+2)
            for i in range(0,self.number_of_fibers):
                self.init_params[2*i]=750 #T1 in ms
                self.lowerbounds[2*i]=200#0
                self.upperbounds[2*i]=4000#np.inf#
                self.init_params[2*i+1]=1.7E-3 #Dpar in mm2/s
                self.lowerbounds[2*i+1]=0.1E-3#0#
                self.upperbounds[2*i+1]=5.0E-3#np.inf




            #additional Johnson noise term neta:
            self.init_params[2*self.number_of_fibers]=0.0 #7.544243 is actual noise mean in asparagus phantom

            #neta:
            self.lowerbounds[2*self.number_of_fibers]=0
            self.upperbounds[2*self.number_of_fibers]=np.inf

            #unknown S0: this depends very much on the acquisition
            #use first signal point (first TI, b=0) and init T1, assume long TR DO change to steady state eq
            self.init_params[2*self.number_of_fibers+1]=np.absolute(self.IR_DWIs[0]/(1-2*np.exp(-1.0*self.TIs[0]/750)))

            #S0:
            self.lowerbounds[2*self.number_of_fibers+1]=0
            self.upperbounds[2*self.number_of_fibers+1]=np.inf





        #we have number of TIs * number of bvalues observations, need to string them all out and add the constants to each observation
        #fastest varying will be diffusion, then TI

        args=np.zeros((self.number_of_TIs*self.number_of_diff_encodes,10+6*self.number_of_fibers))
        #args[:,0]=bvals
        #args[:,1]=observations
        #args[:,2]=TR
        #args[:,3]=not used right now, could use for fixed neta input
        #args[:,4]=TIs
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
        #args[:,15+6*(self.number_of_fibers-1)]=TE

        if (just_b0): #we are just going to repeat the b=0 images
            for i in range(0,self.number_of_TIs):
                #bvals are all zero so leave that as is (init to zero)
                #TIs:
                args[i*self.number_of_diff_encodes:(i+1)*self.number_of_diff_encodes,4]=np.ones(self.number_of_diff_encodes)*self.TIs[i]

        else:

            for i in range(0,self.number_of_TIs):
                #bvals:
                for j in range(0,self.number_of_diff_encodes):
                    #bvals:
                    args[i*self.number_of_diff_encodes+j,0]=self.grad_table.bvals[j]

                    #TIs:
                    args[i*self.number_of_diff_encodes+j,4]=self.TIs[i]




        #constants:
        args[:,5]=self.number_of_fibers*np.ones(self.number_of_TIs*self.number_of_diff_encodes)
        args[:,6]=self.vic*np.ones(self.number_of_TIs*self.number_of_diff_encodes)
        args[:,13+6*(self.number_of_fibers-1)]=self.viso*np.ones(self.number_of_TIs*self.number_of_diff_encodes)
        args[:,14+6*(self.number_of_fibers-1)]=self.B1*np.ones(self.number_of_TIs*self.number_of_diff_encodes)
        for i in range(0,self.number_of_fibers):
            args[:,7+6*i]=self.AFDs[i]*np.ones(self.number_of_TIs*self.number_of_diff_encodes)




        #string out the observations with diffusion fastest varying
        counter=0
        for i in range(0,self.number_of_diff_encodes):
            for j in range(0,self.number_of_TIs):
                if (self.IR_DWIs[i]<9E-9):
                    print "IR DWI image value 0"
                    return None
                if (just_b0): #we will repeat (not necessary but allows the rest of this code to be used):
                    args[j*self.number_of_diff_encodes+i,1]=self.IR_DWIs[j]


                else:

                    args[j*self.number_of_diff_encodes+i,1]=self.IR_DWIs[counter]

                counter+=1



        #the nominal TR and TE:
        args[:,2]=self.TR*np.ones(self.number_of_TIs*self.number_of_diff_encodes)
        args[:,15+6*(self.number_of_fibers-1)]=self.TE*np.ones(self.number_of_TIs*self.number_of_diff_encodes)

        #create the g-vector for each fiber in coord system aligned along that fiber:
        for k in range(0,self.number_of_fibers):



            f=np.zeros(3)

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




        #to weight b=0 more (instead of actually acquiring more), repeat N more times:

        N=np.int(0.1*(self.number_of_diff_encodes-1))-1
        newargs=np.zeros((self.number_of_TIs*self.number_of_diff_encodes+self.number_of_TIs*N,10+6*self.number_of_fibers))
        for i in range(self.number_of_TIs*self.number_of_diff_encodes):
            newargs[i,:]=args[i,:]
        for j in range(self.number_of_TIs):
            for i in range(N):

                newargs[self.number_of_TIs*self.number_of_diff_encodes+j*N+i,:]=args[j*self.number_of_TIs,:]      #the b=0 for that TI



        if (simulate):#simulate data for these input fibers and AFDs:
            #NOTE this is just for Dpar eq case
            #default params, except potentially different T1s and no neta
            new_params=np.copy(self.init_params)
            if set_Dpar_equal:
                new_params[self.number_of_fibers+1]=0.0 #neta
            elif just_b0:
                new_params[1]=0.0 #neta


            sim_data=SignedIRDiffEqn(new_params,newargs)

            #this will just be the magnitude difference:
            #sim_data=newargs[:,1]+IRDiffEqn(new_params,newargs) # this is the equation (sim) result
            random.seed()
            #real only: we don't want to do this
            newargs[:,1]=np.absolute(sim_data+[random.gauss(0,sim_noise_level) for i in range(len(sim_data))])

            #add noise on two channels: we do want to do this

            #newargs[:,1]=np.sqrt(np.square(sim_data+[random.gauss(0,sim_noise_level) for i in range(len(sim_data))])+np.square([random.gauss(0,sim_noise_level) for i in range(len(sim_data))]))
               
        #fit the equation: there are a lot of options here; user can modify this call
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
            #ax.scatter(range(self.number_of_TIs*self.number_of_diff_encodes),args[:,1] , s=2, alpha=0.4)

            #normal 30 dir acq
            #just the b=0, z y x grad orientations
            if (True):


                plotdata=np.zeros([self.number_of_TIs,4])
                thisDWI=[0, 21, 5, 30]
                for i in range(self.number_of_TIs):
                    for j in range(4):
                        plotdata[i,j]=newargs[i*self.number_of_diff_encodes+thisDWI[j],1]
                ax.plot(self.TIs,plotdata[:,0], 'k--')
                ax.plot(self.TIs,plotdata[:,1], 'b--')
                ax.plot(self.TIs,plotdata[:,2], 'g--')
                ax.plot(self.TIs,plotdata[:,3], 'r--')

                ax.set_title('All TIs, b=0 (black), ~z (blue), ~y (green), ~x (red) gradient orientations', fontsize=18)

                #now set to predicted signal:

                pred_sig_res=IRDiffEqn(res_lsq.x,newargs)



                for i in range(self.number_of_TIs):
                    for j in range(4):
                        plotdata[i,j]=newargs[i*self.number_of_diff_encodes+thisDWI[j],1]+pred_sig_res[i*self.number_of_diff_encodes+thisDWI[j]]
                ax.plot(self.TIs,plotdata[:,0], 'k-')
                ax.plot(self.TIs,plotdata[:,1], 'b-')
                ax.plot(self.TIs,plotdata[:,2], 'g-')
                ax.plot(self.TIs,plotdata[:,3], 'r-')

                #these vary a bit with acquisition
                #textstr='dashed: data\nsolid: fit\n\nactual directions:\n~x=[-0.97958797216415, 0.17135678231716, 0.10509198158979]\n~y=[0.20307792723178, 0.94054549932479, -0.27227476239204]\n~z=[-0.20379328727722, 0.17156073451042, 0.96386468410492]'
                #ax.text(0.9,250,textstr)

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
    def SetInputData(self,fiber_dirs,AFDs,IR_DWIs,TIs,grad_table,vic,TR,TE,vox0,vox1,vox2,sag,Dpareq,viso,B1):



        self.AFDs = AFDs

        #each fiber dir has a vector, fiber_dirs is 2D
        self.fiber_dirs = fiber_dirs

        #number_of_fibers=size of AFD array
        self.number_of_fibers=len(self.AFDs)

        self.IR_DWIs = IR_DWIs

        self.TIs = TIs
        print("TIs: ")
        print TIs

        self.number_of_TIs = len(self.TIs)

        self.grad_table = grad_table

        self.number_of_diff_encodes = len(self.grad_table.bvals)/self.number_of_TIs

        self.vic=vic
        self.viso=viso
        self.B1=B1

        self.TR = TR
        self.TE = TE

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




def IRDiffEqn(params,*args): #equation for residuals; params is vector of the unknowns




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
                    if (not set_noIE):
                        term2=SteadyStateT1Recov(params[numfibers+3], B1, TIs[h], TR, TE, params[i])
                    else:
                        term2=SteadyStateT1RecovnoIE(B1, TIs[h], TR, TE, params[i])

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
            if (not set_noIE):
                eq[h]=(1-viso)*eq[h]+viso*(np.exp(-bvals[h]*CSF_D)*SteadyStateT1Recov(params[numfibers+3], B1, TIs[h], TR, TE, CSF_T1))
            else:
                eq[h] = (1 - viso) * eq[h] + viso * (np.exp(-bvals[h] * CSF_D) * SteadyStateT1RecovnoIE(B1, TIs[h], TR, TE, CSF_T1))

        else: #just_b0
            if (not set_noIE):
                term2=SteadyStateT1Recov(params[3], B1, TIs[h], TR, TE, params[0])
            else:
                term2=SteadyStateT1RecovnoIE(B1, TIs[h], TR, TE, params[0])

            #GE ver
            #term2=1-2*np.exp(-1.0*TIs[h]/params[0])+np.exp(-1.0*TR/params[0])

            eq[h]+=term2




        #take magnitude, mult by S0, and add Johnson noise term neta:
        #params[2*numfibers+1] is S0
        #params[2*numfibers] is noise term neta, currently added for ALL images


        if (set_Dpar_equal):
            sig[h]=sqrt((params[numfibers+2]*eq[h])**2+params[numfibers+1]**2)
        elif (just_b0):
            sig[h]=sqrt((params[2]*eq[h])**2+params[1]**2)
        else:  #Dpar not equal, not implemented
            sig[h]=sqrt((params[2*numfibers+1]*eq[h])**2+params[2*numfibers]**2)
        out[h]=sig[h]-obs[h]

    #if (numfibers>1):
        #print("Dpar: %f %f; T1: %f %f; SOS residuals: %f" % (params[1], params[3], params[0], params[2], np.sum(np.square(out))))


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

                    if (not set_noIE):
                        term2=SteadyStateT1Recov(params[numfibers+3], B1, TIs[h], TR, TE, params[i])
                    else:
                        term2=SteadyStateT1RecovnoIE(B1, TIs[h], TR, TE, params[i])

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
            if (not set_noIE):
                eq[h]=(1-viso)*eq[h]+viso*(np.exp(-bvals[h]*CSF_D)*SteadyStateT1Recov(params[numfibers+3], B1, TIs[h], TR, TE, CSF_T1))
            else:
                eq[h]=(1-viso)*eq[h]+viso*(np.exp(-bvals[h]*CSF_D)*SteadyStateT1RecovnoIE(B1, TIs[h], TR, TE, CSF_T1))



        else: #just_b0

            term2=SteadyStateT1Recov(params[3], B1, TIs[h], TR, TE, params[0])

            #GE ver
            #term2=1-2*np.exp(-1.0*TIs[h]/params[0])+np.exp(-1.0*TR/params[0])

            eq[h]+=term2




        #take magnitude, mult by S0, and add Johnson noise term neta:
        #params[2*numfibers+1] is S0
        #params[2*numfibers] is noise term neta, currently added for ALL images


        if (set_Dpar_equal):
            sig[h]=params[numfibers+2]*eq[h]
        elif (just_b0):
            sig[h]=params[2]*eq[h]
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

def SteadyStateT1RecovnoIE(B1, TI, TR, TE, T1):
    #RF pulses
    theta1=1.0*pi
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
        if  iszero(v[3]):
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
    