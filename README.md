# fibermyelin

Estimating orientation-specific T1 values in combined inversion recovery and diffusion (IR-diff) data.

## Pre-processing
Locate all required inputs
`irdiff_preprocess.pl`

## Average tensor estimation and un-shuffling
`irdiff_pipeline.py`
 - Compute an average tensor from HARDI data
 - Un-shuffle the IR-diff data, such that each slice is ordered by TI

## T1 estimation
`fibermyelin_pipeline.py`

<!--
Operation of FiberT1Solver with afd-weighted fiber T1 averaging in a voxel

## To do before running FiberT1Solver:
- Set true_afd_thresh in FiberT1Solver.py (line 60) and in fibermyelin_pipeline.py (line 62). This threshold is equivalent to     the afdthresh we have always been using (i.e. let true_afd_thresh = 0.2).
- Set 'on_true_afd_thresh = True' line 70 in fibermyelin_pipelin.py to activate afd-weighted fiber T1 averaging.
- Set 'afdthresh 0.0' in the fibermyelin_pipeline command. Ensures that a T1 time is calculated for each fiber in each voxel. 

## Fiber-specific T1 calculated as follows:
- Fiber-specific T1 of EACH fiber (even sub-afdthresh ones) is calculated using the SE-IR-DWI eqn in voxels with ANY fiber_afd > true_afd_thresh (happens in FiberT1Solver.py). Calculated T1 times of ALL sub-afdthresh fibers in a voxel are then replaced with the afd-weighted average T1 of the super-afdthresh fibers in that voxel (happens right before T1 times are saved, ~line 423 in fibermyelin_pipeline.py)
- In voxels with ONLY sub-afdthresh fibers, a single T1 time is calculated using the SE-IR equation (no diffusion term, same eqn. as is used for the just_b0 T1 calculation). This T1 time is then assigned to all fibers in the voxel (again before T1 times are saved, ~line 423 in fibermyelin_pipeline.py). 
-->
