
---

# Neural Mechanisms of Resource Allocation in Working Memory

This repository contains data and code reported in:

> **Li, H. H., Sprague, T. C., Yoo, A. H., Ma, W. J., & Curtis, C. E. (2025). Neural mechanisms of resource allocation in working memory. Science Advances.**

If you use any of the code or data in this repository, include the above reference.

---

## I. Data

### Behavioral Data
- **Folders:** `data/behav_1item` and `data/behav_2item`
- **Description:**  
  Files in these folders contain trial-by-trial information on condition, target location, and behavioral reports.
- **Filenames:**
  Files are named by the subject ID (from S1 to S11)
- **Variables for the 1-item Experiment:**
  - `wm_ang`:  
    - Location of the target (in degree polar angle)
  - `behEst`:  
    - Behavioral reports (in polar angle)
- **Variables for the 2-item Experiment:**
  - `targ_angs`:  
    - First column: Location of the target (in degree polar angle)
    - Second column: Location of the non-target (in degree polar angle)
  - `conditions`:  
    - 1 = valid condition (precue matched the target)  
    - 2 = invalid condition
  - `behEst`:  
    - Behavioral reports (in polar angle)
  - `sacRT`:  
    - Saccade reaction time (in second)
  - `goodtrial`:  
    - Index for trials that passed exclusion criteria (1=good trials; 0=excluded trials)
  - `sep_cue_angs`:  
    - Angle (degree) of the separator in the precue indicating how the aperture is divided
  - `pri_cue_angs`:  
    - Angle (degree) of the precue indicating the prioritized item


### fMRI Data
- **Folders:** `data/trialData_1item` and `data/trialData_2item`
- **Description:**  
  Each `.mat` file contains voxel activity pattern (#trial X #voxel) for each subject, session, and ROI. Data from the 1-item experiment is only used for voxel selection.
- **Filenames:**
  Files are named as subjectID_sessionID_ROI_surf_trialData.mat
- **Variables for the 1-item Experiment in `data/trialData_1item`:**
  - `TR`: TR = 0.75 second
  - `c_map`:  
    - Location of the target (in degree polar angle)
  - `dt_mapz`:  
    - Voxel activity patterns (#trial X #voxel)
- **Variables for the 2-item Experiment in `data/trialData_2item`:**
  - `TR`: TR = 0.75 second
  - `c_all`:  
    - First column: Location of the target (in degree polar angle)
    - Second column: Location of the non-target (in degree polar angle)
    - Third column: conditions (1 = valid condition; 2 = invalid condition)
  - `dt_allz`:  
    - Voxel activity patterns (#trial X #voxel)

---

## II. Outputs of Analysis

- **Folder:** `mdata/decoded`
- **Description:**  
  Contains decoding results for each subject and ROI. Each `.mat` file includes
- **Filenames:**
  Files are named as subjectID_sessionID_ROI_decoded_750vox.mat
- **Variables:**  
  - `est`: Decoded locations  
  - `unc`: Decoded uncertainty  
  In both variables, the first column represents the target and the second column represents the distractor (non-target).
  - `p`: Hyperparameters used in decoding
  - `vox_mask`: index for selected voxels
  - `nvox`: number of voxels for decoding
  - `sess_2item`: sessions decoded
  
- **Folder:** `mdata/VPmodel`
- **Description:**  
  Contains fitting results of a priority-dependent variable precision model to the behavioral data (used for Supplementary Fig.2).
- **Variables:**  
  - `berr_mean`: mean absolute behavioral error (valude, invalid)
  - `berr_std`: Sdv of behavioral error (valude, invalid)
  - `bestPar`: best-fit parameters
  - `nLLVec`: negative log-likehood
  - `parMat`: parameters under different intial points

---

## III. Code

- **Main Script:**  
  Run `wmPriority_wrapper.m` to decode two items for each trial from the voxel activity of a particular subject and ROI.

  `wmPriority_genModelDecode_2item.m` is called by the wrapper to read and sort the data for decoding with leave-one-run-out cross validation procedure.
  
- **Decoding Algorithms:**  
  The algorithms and code for conducting Bayesian decoding of 2 items from neural (fMRI) activity are located in the folder `Decode_2item`.

- **Summary:**  
  Run `plot_summary.m` to plot summary of the results shown in the paper. The summary results are saved in summary.mat

- **Others:**
  cat_struct.m, cpsFigure.m, plot_errorbar.m, plot_label.m are miscellaneous functions for plotting or sorting the data

---
