# wmPriority

---

# Neural Mechanisms of Resource Allocation in Working Memory

This repository contains data and code reported in:

> **Li, H. H., Sprague, T. C., Yoo, A. H., Ma, W. J., & Curtis, C. E. (2025). Neural mechanisms of resource allocation in working memory. Science Advances.**

---

## I. Data

### Behavioral Data
- **Folders:**  
  - `data/behav_1item`
  - `data/behav_2item`
- **Description:**  
  These folders contain trial-by-trial information including condition, target location, and behavioral reports.
- **Essential Variables for the 2-item Experiment:**
  - `psy.targ_angs`  
    * First column: polar angle of the (saccade) target  
    * Second column: polar angle of the distractor (non-target)
  - `psy.conditions`  
    * 1 = valid condition (precue matched the target)  
    * 2 = invalid condition
  - `psy.behEst`  
    * Behavioral reports (in polar angle, readout from the eye-tracking data)

### fMRI Data
- **Folders:**  
  - `data/trialData_1item`
  - `data/trialData_2item`
- **Description:**  
  Each `.mat` file contains voxel activity pattern (#trial X #voxel) for each subject, session, and ROI. Data from the 1-item experiment is used only for voxel selection.

### Outputs of Analysis
- **Folder:** `mdata/decoded`
  - **Contents:**  
    Decoding results for each subject and ROI. Each `.mat` file contains:
    - `lf`: Decoded 2-dimensional likelihood surface
    - `est`: Decoded locations (first column: target, second column: distractor)
    - `unc`: Decoded uncertainty (first column: target, second column: distractor)
    - `p`: Hyperparameters used in decoding (includes conditions, target polar angle, and non-target polar angle)
- **Folder:** `mdata/VPmodel`
  - **Contents:**  
    Fitting results of a priority-dependent variable precision model to the behavioral data (used for Supplementary Fig.2)

---

## II. Code

- **Main Script:**  
  Run `wmPriority_wrapper.m` to decode two items from the voxel activity of a particular subject and ROI.
- **Decoding Algorithms:**  
  Located in the folder `Decode_2item`. The estimation of the noise covariance matrix adapts methods and code developed by [the original authors or source, if applicable].

---
