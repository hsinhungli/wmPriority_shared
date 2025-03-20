This repository contains data and code reported in 

Li, H. H., Sprague, T. C., Yoo, A. H., Ma, W. J., & Curtis, C. E. (2025). Neural mechanisms of resource allocation in working memory. Science Advances.

I. Data:
Behavioral data are in folders data/behav_1item and data/behav_2item.
Files in these folders contain trial-by-trial information of the condition, target location and behavioral reports. Below are some variables essential to the 2-item experiment:
psy.targ_angs : the first column is the polar angle of the (saccade) target, and the second column is the distractor (non-target)
psy.conditions: 1 = valid condition (precue matched the target); 2 = invalid condition
psy.behEst    : Behavioral reports (in polar angle, readout from the eye-tracking data)

fMRI voxel activity used for decoding is in folders data/trialData_1item and data/trialData_2item. 
Each mat file contains voxel activity pattern (#trial X #voxel) for each subject, session, and ROI. 
Data of the 1-item experiment is only used for voxel selection. 

II. Outputs of analysis:
mdata/decoded contains decoding results for each subject and ROI. Decoded locations and decoded uncertainty are stored in variables "est" and "unc" respectively, in which the first column represents the target and the second column represents the distractor (non-target).
Variable "p" contains some hyperparameters used in decoding. For convenience, the conditions (valid, invalid condition), the polar angle of the target and the non-target are also stored in "p".
mdata/VPmodel contains fitting results of a priority-dependent variable precision model to the behavioral data (used for Supplementary Fig.2)

III. Code
Run wmPriority_wrapper.m to decode two items the voxel activity of a particular subject and ROI.
The decoding algorithms are in the folder Decode_2item. The estimation of the noise covariance matrix adapted the TAFKAP method developed by Van Bergen, JFM Jehee - BioRxiv, 2021.
Run plot_summary.m to plot summary of the results shown in the paper.
