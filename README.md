
---

# Neural Mechanisms of Resource Allocation in Working Memory

This repository contains data and code reported in:

> **Li, H. H., Sprague, T. C., Yoo, A. H., Ma, W. J., & Curtis, C. E. (2025). Neural mechanisms of resource allocation in working memory. Science Advances.**

---

## I. Data

### Behavioral Data
- **Folders:** `data/behav_1item` and `data/behav_2item`
- **Description:**  
  Files in these folders contain trial-by-trial information on condition, target location, and behavioral reports.
- **Essential Variables for the 2-item Experiment:**
  - `psy.targ_angs`:  
    - First column: polar angle of the (saccade) target  
    - Second column: polar angle of the distractor (non-target)
  - `psy.conditions`:  
    - 1 = valid condition (precue matched the target)  
    - 2 = invalid condition
  - `psy.behEst`:  
    - Behavioral reports (in polar angle, readout from the eye-tracking data)

### fMRI Data
- **Folders:** `data/trialData_1item` and `data/trialData_2item`
- **Description:**  
  Each `.mat` file contains voxel activity pattern (#trial X #voxel) for each subject, session, and ROI. Data from the 1-item experiment is only used for voxel selection.

---

## II. Outputs of Analysis

- **Folder:** `mdata/decoded`
  - **Description:**  
    Contains decoding results for each subject and ROI. Each `.mat` file includes:
    - **Variable `lf`:** Decoded 2-dimensional likelihood surface.
    - **Variables `est` and `unc`:**  
      - `est`: Decoded locations  
      - `unc`: Decoded uncertainty  
      In both variables, the first column represents the target and the second column represents the distractor (non-target).
    - **Variable `p`:**  
      Contains hyperparameters used in decoding, including conditions (valid, invalid), the polar angle of the target, and the non-target.
  
- **Folder:** `mdata/VPmodel`
  - **Description:**  
    Contains fitting results of a priority-dependent variable precision model to the behavioral data (used for Supplementary Fig.2).

---

## III. Code

- **Main Script:**  
  Run `wmPriority_wrapper.m` to decode two items from the voxel activity of a particular subject and ROI.
  
- **Decoding Algorithms:**  
  Located in the folder `Decode_2item`. The estimation of the noise covariance matrix adapts the TAFKAP method developed by Van Bergen & Jehee - *BioRxiv, 2021*.

---
