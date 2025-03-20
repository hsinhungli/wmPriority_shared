function [est1, est2, unc1, unc2] = Readout(lf_surf, psy, p)

%This function read out two decoded locations based on the 2D likelihood surface. 
%lf_surf: should be three dimensional: n_angs X n_angs X ntrial
%During decoding, for each likelihood surface, lf_surf(:,:,tt), we define the horizontal axis as the
%high-priority item, and the vertical axis as the low-priority item (weighted by a gain factor)
%psy: is trial-by-trial task-related information. During decoding we
%utilize the knowledge of cue direction, and which half of the aperture is
%cued as high-priority

%n_angs = 360;
%angsd = linspace(360/n_angs,360,n_angs); %in degree
%angs = 2*pi*angsd/360; %in radian
%[ang1d,ang2d] = meshgrid(angsd,angsd);

angsd = p.angsd;
angs = p.angs;
[ang1d,ang2d] = meshgrid(angsd,angsd);
ntrial = size(lf_surf, 3);
angle_diff = @(a,b) mod((a-b) + 180, 360)-180;

%We need to use this later. We assume that the decoder knows which side of the aperture is cued as high-priority
high_loc = psy.targ_angs(sub2ind([ntrial 2],(1:ntrial)',psy.conditions));

est1 = nan(ntrial,1);
unc1 = nan(ntrial,1);
est2 = nan(ntrial,1);
unc2 = nan(ntrial,1);

pri_cue_angs = mod(mod(psy.pri_cue_angs, 360), 180); %this is the direction of the precue
for tt = 1:ntrial
    thislf = lf_surf(:,:,tt);
    thiscue = pri_cue_angs(tt);
    
    %implementing task-constraint like a prior: horizontal dimension is stimulus CW to the separator
    prior = angle_diff(ang1d, thiscue) < 0 & angle_diff(ang2d, thiscue) >= 0; 
    
    % In lf_surf, the first (vertical) dimension is the low-priority. The horizontal dimension is high-priority.
    % We assume that the decoder knows which side of the aperture is cued as high-priority. 
    % Thus, flip the x- and y-axis of lf_surf if the high-priority item is NOT CW to the separator
    flipInd = angle_diff(high_loc(tt), thiscue) > 0;
    if flipInd == 1
        thislf = thislf';
    end
    posterior = thislf .* prior;
    
    lf_1 = squeeze(sum(posterior,1));
    lf_2 = squeeze(sum(posterior,2));
    
    lf_1 = lf_1 ./ sum(lf_1);
    lf_2 = lf_2 ./ sum(lf_2);
    
    pop_vec = lf_1*exp(1i*angs');
    est1(tt) = mod(angle(pop_vec)/pi*180, 360); %Stimulus estimate (likelihood/posterior means)
    unc1(tt) = sqrt(-2*log(abs(pop_vec)))/pi*180;
    
    pop_vec = lf_2'*exp(1i*angs');
    est2(tt) = mod(angle(pop_vec)/pi*180, 360); %Stimulus estimate (likelihood/posterior means)
    unc2(tt) = sqrt(-2*log(abs(pop_vec)))/pi*180;
    
    %Save 1-D posterior if needed
    %lf_1Mat(tt,:) = lf_1;
    %lf_2Mat(tt,:) = lf_2;
end

target_label = angle_diff(psy.targ_angs, pri_cue_angs) < 0;
flipInd = target_label(:,1) == 0; %flip if the target is NOT CW to the separator
est1_orig = est1(flipInd);
est2_orig = est2(flipInd);
unc1_orig = unc1(flipInd);
unc2_orig = unc2(flipInd);

est1(flipInd) = est2_orig;
est2(flipInd) = est1_orig;
unc1(flipInd) = unc2_orig;
unc2(flipInd) = unc1_orig;

% lf_1Mat_orig = lf_1Mat(flipInd,:);
% lf_2Mat_orig = lf_2Mat(flipInd,:);
% lf_1Mat(flipInd,:) = lf_2Mat_orig;
% lf_2Mat(flipInd,:) = lf_1Mat_orig;

return

end