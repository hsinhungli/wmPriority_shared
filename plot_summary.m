subjPool = {'S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11'};
ROIs = {'V1','V2','V3','V3AB','IPS0','IPS1','IPS2','IPS3','sPCS','iPCS'};
nsub = length(subjPool);
nROI = length(ROIs);
mld = @(A,B) ([A ones(size(A,1),1)])\B;
load('summary.mat');
%% Figure 4: Decoding error
cpsFigure(1,.7);
data = sq(mdata.nerr_mean(1,:,:));
myerrorbar((1:nROI)-.1, mean(data,2), std(data,[],2)/sqrt(nsub),0,[.2 .2 .2]);
data = sq(mdata.nerr_mean(2,:,:));
myerrorbar((1:nROI)+.1, mean(data,2), std(data,[],2)/sqrt(nsub),0,[.7 .7 .7]);
xlim([.5 nROI+.5])
ylim([20 60]);
set(gca, 'XTick', 1:nROI, 'XTickLabel', ROIs, 'YTick', 20:10:60, 'FontSize',12);
box off;
ylabel('Decoding Error');

%% Figure 5: Gains for the low-priority items
cpsFigure(1,.5);
data = mdata.est_w; 
thisx = 1:nROI;
thisy = mean(data,2);
sem = std(data,[],2) / sqrt(nsub);
for vv = 1:nROI
    ss = scatter(ones(1,nsub)*vv, data(vv,:),45,'filled'); hold on;
    ss.MarkerEdgeColor = [1 1 1];
    ss.MarkerFaceColor = [.8 .8 .8];
end
myerrorbar(thisx, thisy, sem, 0,'k',[]); hold on;
ff = fplot(@(x) x*0+1, [0 nROI+1],'--');
ff.Color = [.5 .5 .5];
axis([.5 nROI+.5 -.1 2.1]);
box off;
set(gca, 'XTick',1:nROI,'XTickLabel',ROIs,'FontSize',14);
ylabel('gain of the low-priority item','FontSize',14);

%% Figure 6A: Decoded uncertainty
cpsFigure(1,.7);
data = sq(mdata.unc(1,:,:));
myerrorbar((1:nROI)-.1, mean(data,2), std(data,[],2)/sqrt(nsub),0,[.2 .2 .2]);
data = sq(mdata.unc(2,:,:));
myerrorbar((1:nROI)+.1, mean(data,2), std(data,[],2)/sqrt(nsub),0,[.7 .7 .7]);
axis([ .5 nROI+.5 10 50]);
set(gca, 'XTick', 1:nROI, 'XTickLabel', ROIs, 'FontSize',12);
box off;
ylabel('Decoded uncertainty (deg)');

%% Figure 6B. Binned Correlation: Neu Pri vs. Beh Pri (berr)
cpsFigure(1,.8);
cmap = parula(6);
pval = nan(nROI, 1);
for vv = 1:nROI
    subplot(2,ceil(nROI/2),vv)
    nbin = mdata.nbin(1);
    xdata = sq(mdata.preAMI_bin(:,vv,:));
    ydata = -sq(mdata.bPri_preAMIbin(:,vv,:));
    
    xerr = xdata - mean(xdata) + mean(xdata(:));
    yerr = ydata - mean(ydata) + mean(ydata(:));
    xerr = reshape(xerr,size(xdata));
    yerr = reshape(yerr,size(ydata));
    
    thiscorr = corr(xerr(:), yerr(:));
    cf_temp = mld(xerr(:), yerr(:));
    fplot(@(x) cf_temp(1).*x+cf_temp(2), [-.8 1],'LineWidth',2, 'Color', [.5 .5 .5 .8]); hold on;
    
    for bb = 1:nbin
        ss = scatter(xerr(bb,:),yerr(bb,:),'filled'); hold on;
        ss.MarkerFaceColor = cmap(bb+1,:);
        ss.MarkerFaceAlpha = .8;
        ss.MarkerEdgeColor = [1 1 1];
    end
    
    title(ROIs{vv});
    text(-1,-6.5,sprintf('Corr = %1.2f',thiscorr),'FontSize',10);
    
    if vv==1
    mylabel('Neural prioritization', 'Behavioral prioritization');
    end
    axis([-1 1 -8 8]);
    axis square;
    box off;
end
tightfig;

%% Figure 6D. Binned Correlation: Neu Pri vs. Beh Pri (sacRT)
cpsFigure(1,.8);
cmap = parula(6);
for vv = 1:nROI
    subplot(2,ceil(nROI/2),vv)
    nbin = 4;
    
    xdata = sq(mdata.nDiff_bin(:,vv,:));
    ydata = sq(mdata.nDiff_sacRT_bin(:,vv,:));
    
    xerr = xdata - mean(xdata) + mean(xdata(:));
    yerr = ydata - mean(ydata) + mean(ydata(:));
    xerr = reshape(xerr,size(xdata));
    yerr = reshape(yerr,size(ydata));
    
    thiscorr = corr(xerr(:), yerr(:));
    cf_temp = mld(xerr(:), yerr(:));
    fplot(@(x) cf_temp(1).*x+cf_temp(2), [-1 max(xdata(:))*1.1],'LineWidth',2, 'Color', [.5 .5 .5 .8]); hold on;
    
    for bb = 1:nbin
        ss = scatter(xerr(bb,:),yerr(bb,:),40,'filled'); hold on;
        ss.MarkerFaceColor = cmap(bb+1,:);
        ss.MarkerFaceAlpha = .8;
        ss.MarkerEdgeColor = [1 1 1];
    end

    title(ROIs{vv});
    text(-.5,.475,sprintf('Corr = %1.2f',thiscorr),'FontSize',10);
    
    if vv==1
    mylabel('Probe prioritization', 'Saccade RT (sec)');
    end
    axis([-1 1 .45 .65]);
    axis square;
    box off;
end

%% Figure S1. Decoding Sdv
cpsFigure(1,.7);
data = sq(mdata.nerr_std(1,:,:));
myerrorbar((1:nROI)-.1, mean(data,2), std(data,[],2)/sqrt(nsub),0,[.2 .2 .2]);
data = sq(mdata.nerr_std(2,:,:));
myerrorbar((1:nROI)+.1, mean(data,2), std(data,[],2)/sqrt(nsub),0,[.7 .7 .7]);
xlim([.5 nROI+.5])
set(gca, 'XTick', 1:nROI, 'XTickLabel', ROIs, 'FontSize',12);
box off;
ylabel('Decoding Error');

%% SupFig2. Across Subject - Neural prioritization vs. weight in VP model
cpsFigure(1.6,.8);
for ss = 1:nsub 
    fn2l = sprintf('mdata/VPModel/%s_bestPar.mat',subjPool{ss});
    thisPar = load(fn2l);
    Jbar(ss) = thisPar.bestPar(1);
    tau(ss) = thisPar.bestPar(2);
    pweight(ss) = thisPar.bestPar(3);
end

preAMI = sq(mean(mdata.preAMI_bin, 1));
for vv = 1:nROI
    subplot(2,ceil(nROI/2),vv);
    thisy = pweight';
    thisx = preAMI(vv,:)';

    cf_temp = mld(thisx(:), thisy(:));
    fplot(@(x) cf_temp(1).*x+cf_temp(2), [min(thisx) max(thisx)],'LineWidth',2, 'Color', [.8 .8 .8 .5]); hold on;
    [thiscorr, ~]=corr(thisx(:), thisy(:));

    title(ROIs{vv});
    ss = scatter(thisx, thisy, 'filled');
    ss.MarkerFaceAlpha = .8;
    ss.MarkerFaceColor = [0 .3 1];
    text(.4,.3,sprintf('Corr = %1.2f',thiscorr),'FontSize',10);
    axis([-.5 1 0.2 1]);
    
    if vv==1
        mylabel('Neural prioritization','Estimated weight',12);
    end
    
    box off; 
end

%% Figure S3. Across Subjects: Neural prioritization vs. Saccade RT Difference
cpsFigure(1,.8);
for vv = 1:nROI
    
    subplot(2,ceil(nROI/2),vv);
    thisy = mdata.sacRT(2,:) - mdata.sacRT(1,:);
    
    nPri = sq(mean(mdata.preAMI_bin,1));
    thisx = nPri(vv,:)';
    cf_temp = mld(thisx(:), thisy(:));
    fplot(@(x) cf_temp(1).*x+cf_temp(2), [min(thisx) max(thisx)],'LineWidth',2, 'Color', [.8 .8 .8 .5]); hold on;
    [thiscorr, pval(vv)]=corr(thisx(:), thisy(:));
    title(ROIs{vv});
    text(.02,-.02,sprintf('Corr = %1.2f',thiscorr),'FontSize',10);
    
    ss = scatter(thisx, thisy, 'filled');
    ss.MarkerFaceAlpha = .8;
    ss.MarkerFaceColor = [0 .3 1];
    axis([-.5 1 -.05 .2]);
    axis square
    box off;
    
    if vv==1
    mylabel('Unc L-H','sacRT L-H',12);
    end
    
end

%% Figure S4. Decoded uncertainty vs. RT - single-trial correlation
cpsFigure(1,.7);
data = sq(mdata.corr_unc_sacRT_pool);
thisx = 1:nROI;
thisy = mean(data,2);
sem = std(data,[],2) / sqrt(nsub);
for vv = 1:nROI
    ss = scatter(ones(1,nsub)*vv, data(vv,:),'filled'); hold on;
    ss.MarkerEdgeColor = [1 1 1];
    ss.MarkerFaceColor = [.8 .8 .8];
end
myerrorbar(thisx, thisy, sem, 0,'k',[]); hold on;
ff = fplot(@(x) 0, [0 nROI+1],'--');
ff.Color = [.5 .5 .5];
axis([.5 nROI+.5 -.4 .4]);
box off;
set(gca,'XTick',1:nROI, 'XTickLabel',ROIs,'YTick',-.3:.1:.3,'FontSize',14);
ylabel('Correlation','FontSize',14);
