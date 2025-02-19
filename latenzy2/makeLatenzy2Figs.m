function figHandles = makeLatenzy2Figs(sLatenzy2,spikeTimes1,eventTimes1,spikeTimes2,eventTimes2,useMaxDur,makePlots)
% make figures for latenzy2 method
%
% history:
%   v0.9 - 19 February 2025
%   - created by Robin Haak

%% prep
%#ok<*AGROW>
latency = sLatenzy2.latency;
peakTimes = sLatenzy2.peakTimes;
peakVals = sLatenzy2.peakVals;
realFrac = sLatenzy2.realFrac;
fracDiff = sLatenzy2.fracDiff;
fracLin = sLatenzy2.fracLin;
realDiff = sLatenzy2.realDiff;
realTime = sLatenzy2.realTime;
meanRealDiff = sLatenzy2.meanRealDiff;
randDiff = sLatenzy2.randDiff;
randTime = sLatenzy2.randTime;
meanRandDiff = sLatenzy2.meanRandDiff;
peakZ = sLatenzy2.peakZ;
latenzyIdx = sLatenzy2.latenzyIdx;

numIters = numel(peakTimes);
useColors = getAbyss(numIters);
lineWidth = 1.5;
markerSize = 60; %50;

%% PLOT
figure;
figHandles = nan(1,6);
if makePlots==1
    %raster plot condition 1
    figHandles(1) = subplot(2,3,1);
    rasterPlot(spikeTimes1,eventTimes1,useMaxDur);
    yLim = ylim;
    hold on
    plot([latency latency], yLim,'color',[0.8627 0.0784 0.2353],'LineStyle','--','LineWidth',lineWidth);
    set(gca,'box','off','TickDir','out');
    ylabel('Trial');
    xlabel('Time from event (s)');
    title('Raster plot condition 1');

    %raster plot condition 2
    figHandles(4) = subplot(2,3,4);
    rasterPlot(spikeTimes2,eventTimes2,useMaxDur);
    yLim = ylim;
    hold on
    plot([latency latency], yLim,'color',[0.8627 0.0784 0.2353],'LineStyle','--','LineWidth',lineWidth);
    set(gca,'box','off','TickDir','out');
    ylabel('Trial');
    xlabel('Time from event (s)');
    title('Raster plot condition 2');
end

%plot cumulative spikes over time
figHandles(2) = subplot(2,3,2); hold on
p = [];
p(1) = plot(realTime{1},realFrac{1,1},'color',[0.8510 0.4510 0.1340],'LineWidth',lineWidth);
p(2) = plot(realTime{1},realFrac{2,1},'color',[0.9255 0.7255 0.5670],'LineWidth',lineWidth);
set(gca,'box','off','TickDir','out');
xlim(useMaxDur);
xlabel('Time from event (s)');
ylabel('Cumulative num. spikes/rep.');
title('Cumulative spikes');
lgd = legend(p,{'1','2'},'Location','southeast','Box','off');
title(lgd, 'Condition');

%plot difference
figHandles(3) = subplot(2,3,3); hold on
p = [];
for iter = 1:numIters
    plot(realTime{iter},fracLin{iter},'color',[0.5 0.5 0.5],'LineWidth',lineWidth);
    p(iter) = plot(realTime{iter},fracDiff{iter},'color',useColors(iter,:),'LineWidth',lineWidth);
    labels{iter} = sprintf('%d', iter);
end
set(gca,'box','off','TickDir','out');
xlim(useMaxDur);
xlabel('Time from event (s)');
ylabel('Deviation (XXXXXXXX)');
title('Cumulative difference');
lgd = legend(p, labels, 'Location', 'best','Box','off');
title(lgd, 'Iteration');

%plot offset from linear baseline
figHandles(5) = subplot(2,3,5); hold on
plot(useMaxDur,[0 0],'color',[0.5 0.5 0.5],'LineWidth',lineWidth,'LineStyle','--');
% p = [];
for iter = 1:numIters
    p(iter) = plot(realTime{iter},realDiff{iter},'color',useColors(iter,:),'LineWidth',lineWidth);
    labels{iter} = sprintf('%d', iter);
end
scatter(peakTimes(~latenzyIdx),peakVals(~latenzyIdx),markerSize,'x','MarkerEdgeColor',[0 0 0],'LineWidth',lineWidth);
scatter(peakTimes(latenzyIdx),peakVals(latenzyIdx),markerSize,'x','MarkerEdgeColor',[0.8627 0.0784 0.2353],'LineWidth',lineWidth);
set(gca,'box','off','TickDir','out');
xlim(useMaxDur);
xlabel('Time from event (s)');
ylabel('Deviation (Δfraction)');
title('Deviation from linear, \color[rgb]{0.8627,0.0784,0.2353}x\color{black} = estimated latency');
% lgd = legend(p, labels, 'Location', 'southeast','Box','off');
% title(lgd, 'Iteration');

%for iteration w/latency, plot real + jitters
figHandles(6) = subplot(2,3,6); hold on
for thisShuffle=1:length(randDiff)
    plot(randTime{thisShuffle,latenzyIdx},(randDiff{thisShuffle,latenzyIdx}-meanRandDiff(thisShuffle,latenzyIdx)) ,'Color',[0.5 0.5 0.5 0.5],'LineWidth',lineWidth);
end
plot(realTime{latenzyIdx},(realDiff{latenzyIdx}-meanRealDiff(latenzyIdx)),'color',useColors(latenzyIdx,:),'LineWidth',lineWidth);
scatter(peakTimes(latenzyIdx),(peakVals(latenzyIdx)-meanRealDiff(latenzyIdx)),markerSize,'x','MarkerEdgeColor',[0.8627 0.0784 0.2353],'LineWidth',lineWidth);
if latenzyIdx(1),xlim(useMaxDur);
else, xlim([useMaxDur(1)  peakTimes(find(latenzyIdx)-1)]);end
set(gca,'box','off','TickDir','out');
xlabel('Time from event (s)');
ylabel('Deviation (Δfraction)');
title(sprintf('Real data + resamplings, mean-subtracted (Z=%.1f)',peakZ(latenzyIdx)));

%add title
sgtitle(sprintf('Estimated latency of the difference = %.4fs', latency), 'FontWeight', 'bold');
end