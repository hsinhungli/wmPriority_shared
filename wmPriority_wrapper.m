addpath(genpath('Decode_2item'));

subjectIndex = 1; %Which subject to decode
roiIndex = 1; %Which ROI

subj = {'S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11'};

%sessions for each subject - main experiment
sess = {{'wmPri1','wmPri2'},{'wmPri1','wmPri2'},{'wmPri1','wmPri2'},{'wmPri1','wmPri2'},{'wmPri1','wmPri2'},...
    {'wmPri1','wmPri2'},{'wmPri1','wmPri2'},{'wmPri1','wmPri2'},{'wmPri1','wmPri2'},{'wmPri1','wmPri2','wmPri3'},...
    {'wmPri1','wmPri2'}};

%sessions for each subject - 1-item experiment for voxel selection
sess_1item = {{'MGSMap1','MGSMap2'},{'MGSMap1','MGSMap2'},{'MGSMap1','MGSMap2'},{'MGSMap1','MGSMap2','MGSMap3'},{'MGSMap1','MGSMap2'},...
    {'MGSMap1','MGSMap2'},{'MGSMap1','MGSMap2'},{'MGSMap1','MGSMap2'},{'MGSMap1','MGSMap2'},{'MGSMap1','MGSMap2','MGSMap3'},...
    {'MGSMap1','MGSMap2'}};

ROIs = {'V1','V2','V3','V3AB','IPS0','IPS1','IPS2','IPS3','iPCS','sPCS'};

savefile = 1; %if set as 1, decoding results will be saved under 'mdata' folder
nvox = 750;

fprintf('\n\n ---------- TAFKAP Starting Subj %s ROI %s nVox %d---------- \n\n', subj{subjectIndex}, ROIs{roiIndex}, nvox);
[lf, est, unc, hypers, p] = wmPriority_genModelDecode_2item(subj{subjectIndex}, sess_1item{subjectIndex}, sess{subjectIndex}, ROIs{roiIndex}, nvox, savefile);

% Outputs: 
% lf: decoded 2D likelihood surface
% est: decoded location (first column is the saccade target ? probed item)
% unc: decoded uncertainty (first column is the saccade target ? probed item)
% hypers: hyperparameters
% p: setting for the model, and some task-related info (p.stimpos is the true locations of the two items)

%% to check decoding results of a subject-ROI: 
% fn2l = 'mdata/decoded/S1_wmPri1wmPri2_V1_decoded_750vox.mat';
% load(fn2l);
% figure; scatter(p.stimpos(:,1), est(:,1)); xlabel('target'); ylabel('decoded location');
% figure; scatter(p.stimpos(:,2), est(:,2)); xlabel('non-target'); ylabel('decoded location');
