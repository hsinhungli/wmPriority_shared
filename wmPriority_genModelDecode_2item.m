function [lf, est, unc, hypers, p] = wmPriority_genModelDecode_2item(subj,sess_1item,sess_2item,ROI,nvox,savefile)

% This function (1) performs voxel selection based on a dataset where the
% subjects performed a 1-item spatial VWM task. (2) Load the voxel activity
% pattern from the 2-item spatial VWM task (the main experiment) and start
% Bayesian decoding with leave-one-run-out cross validation using the
% Decode2item function.

%Inputs: 
%subj: subject ID;
%sess_1item: session names for the 1-item dataset
%sess_2item: session names for the 2-item dataset
%ROI: ROI to decode
%nvox: number of voxel to select
%savefile: flag 1 or 0, to save file or not

if nargin < 6
    savefile = 0;
end
addpath(genpath('./Decode_2item'));

%% Voxel selection load 1-item dataset
fn_psy = sprintf('data/behav_1item/%s_log.mat',subj);
fprintf('BEHAVIOR: %s...\n',fn_psy);
psy = load(fn_psy); %load behavioral data saved in structure 'psy'

all_data = [];
for sess_idx = 1:length(sess_1item)
    
    fn = sprintf('data/trialData_1item/%s_%s_%s_surf_trialData.mat',subj,sess_1item{sess_idx},ROI);
    fprintf('loading %s...\n',fn);
    thisdata = load(fn);
    
    data = [];
    data.dt_z = thisdata.dt_mapz;
    all_data = cat_struct(all_data, data);
    clear data;
end

wm_ang_1item = psy.wm_ang;
data_1item = all_data.dt_z; %% trial X nVoxel
clear all_data psy;

%% load 2-item dataset
%load subj behavioral data
fn_psy = sprintf('data/behav_2item/%s_log.mat',subj);
fprintf('BEHAVIOR: %s...\n',fn_psy);
psy = load(fn_psy); %load behavioral data saved in structure 'psy'

startidx = 1;
all_data = [];
for sess_idx = 1:length(sess_2item)
    
    fn = sprintf('data/trialData_2item/%s_%s_%s_surf_trialData.mat',subj,sess_2item{sess_idx},ROI);
    fprintf('loading %s...\n',fn);
    data = load(fn);
    assert(length(psy.runNumber(psy.sessionNumber==sess_idx))==size(data.dt_allz,1)); %making sure the trial length of behavioral data matches the fMRI data
    data.run_ind = psy.runNumber(psy.sessionNumber==sess_idx);
    all_data = cat_struct(all_data, data);
    all_data.run_ind(startidx:end) = all_data.run_ind(startidx:end) + 100*sess_idx;
    startidx = size(all_data.run_ind,1)+1;
    clear data;
end
all_data.sess_ind = floor(all_data.run_ind/100);
[~,~,all_data.run_ind] = unique(all_data.run_ind);
mydata_trn = all_data.dt_allz; %ntrial X nVoxel

%% select voxels
mystd = std(data_1item,[],1);

% below: adapted from MGSMap_channelRespAmp functions
if nvox < 1 %use all the voxels; only check in case there are bad/dead voxels
    vox_mask = mystd~=0 & ~isnan(mystd); %
    
else %otherwise, we're using the top N voxels sorted by quadrant-wise F-score
    
    vox_mask = zeros(1, size(mydata_trn,2)); %Matrix to save good voxel
    goodvox = mystd~=0 & ~isnan(mystd);
    
    %voxel selection using ANOVA to classify the quadrant of the target
    allF = nan(size(data_1item,2),1);
    allp = nan(size(data_1item,2),1);
    [~,~,thisG] = unique(wm_ang_1item(:));
    
    for voxidx = 1:size(data_1item,2)
        thisX = data_1item(:,voxidx);
        [p,tab,~] = anova1(thisX,thisG,'off');
        allF(voxidx) = tab{2,5};
        allp(voxidx) = p;
        clear thisX p tab stats;
    end
    
    % get rid of NaN's, which can in principle happen when
    % zero std dev, etc.
    allF(isnan(allF)) = -Inf; allF(~goodvox) = -Inf;
    [~, fIdx] = sort(allF,'descend');
    
    rmIdx = ismember(fIdx, find(~goodvox)); %remove bad voxels
    fIdx(rmIdx) = [];
    
    if nvox > length(fIdx)
        these_vox = find(goodvox); %good voxels less than nvoxels wanted, include all then.
    else
        these_vox = fIdx(1:nvox);
    end
    vox_mask(1, these_vox) = 1;
end
fprintf('data Matrix: %s X %s ---\n',num2str(size(mydata_trn,1)),num2str(mean(sum(vox_mask,2))));

%% Setting for training the model
p.stimpos = psy.targ_angs; %target location; first column is the sacade target (probed item)
p.condition = psy.conditions;
p.goodIdx = psy.goodtrial(end-length(all_data.run_ind)+1:end);
p.nchan = 8; %number of basis used in modeling voxel response
p.runNs = all_data.run_ind;
p.wlevel = linspace(0,2,200);
p.ntrial = size(mydata_trn,1);

p.n_angs = 360;
p.angsd = linspace(360/p.n_angs,360,p.n_angs);
p.angs = 2*pi*p.angsd/360;
[p.ang1, p.ang2]  = meshgrid(p.angs, p.angs);
%%
nrun = max(all_data.run_ind);
lf   = nan(360,360,p.ntrial);
hypers = nan(nrun,3); %save hyperparameters

for testrun_idx = 1:nrun % parfor testrun_idx = 1:nrun (implement parallel when running multiple cpus.)
    thisp = p;
    thisp.test_trials = all_data.run_ind == testrun_idx;
    thisp.train_trials = all_data.run_ind ~= testrun_idx;
    
    fprintf(1,'\n\n ---------- Start testrun Idx %s Subj %s ROI %s---------- \n\n', num2str(testrun_idx), subj, ROI);
    [this_lf, this_hypers] = Decode2item(mydata_trn(:,vox_mask==1), thisp);
    fprintf(1,'\n\n ---------- Done testrun Idx %s Subj %s ROI %s ----------\n\n', num2str(testrun_idx), subj, ROI);
    
    temp_lf{testrun_idx} = this_lf;
    temp_hyper{testrun_idx} = this_hypers;
end

%Put decoded results into the shape we want
for testrun_idx = 1:nrun
    lf(:,:,all_data.run_ind == testrun_idx) = temp_lf{testrun_idx};
    hypers(testrun_idx,:) = temp_hyper{testrun_idx};
end

[est1, est2, unc1, unc2] = Readout(lf, psy, p);
est = [est1 est2];
unc = [unc1 unc2];

if savefile == 1
    %Save the decoded likelihood surface
    fn2s = sprintf('mdata/%s_%s_%s_decoded_%ivox.mat',subj,horzcat(sess_2item{:}),ROI,nvox);
    save(fn2s,'est','unc','p','sess_2item','vox_mask','hypers','nvox');
    fprintf(1,'\n ----------Saved filename %s ----------\n',fn2s);
end
end
