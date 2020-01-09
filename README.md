# fibermyelin

# Operation of FiberT1Solver with fiber T1 averaging in the voxel (branch SE-no-IE_voxel_T1_averaging)

## To do before running FiberT1Solver:
- Set true_afd_thresh in FiberT1Solver.py (line 60) and in fibermyelin_pipeline.py (line 62). This threshold is equivalent to     the afdthresh we have always been using.
- Set on_true_afd_thresh on line 70 in fibermyelin_pipelin.py to activate true_afd_thresh-dependent fiber T1 averaging.
- Always have afdthresh = 0.0 in the fibermyelin_pipeline command. Ensures that a T1 time is calculated for each fiber in each voxel. T1 averaging for sub_afd_thresh fibers is calculated later in the script as described below. 

## FiberT1Solver now calculates fiber-specific T1 as follows:
- Fiber-specific T1 of EACH fiber is calculated in voxels with ANY fiber_afd > true_afd_thresh (happens in FiberT1Solver.py). T1 of ALL sub_afd_thresh fibers in a voxel are then set to the afd-weighted average T1 of the super_afd_thresh fibers in that voxel (happens right before T1 times are saved, ~line 423 in fibermyelin_pipeline.py)
- In voxels with ONLY sub_afd_thresh fibers, a single T1 time is calculated using the SE-IR equation (no diffusion term, same equation as is used for the just_b0 T1 calculation). This T1 time is then assigned to all fibers in the voxel. 
