function figHandles = makeLatenzyFigs(sLatenzy,spikeTimes,eventTimes,useMaxDur,makePlots)
% make figures for latenzy method
% Robin Haak, 2024
%
% history:
% 9 August 2024
% - created by Robin Haak
% 16 September 2024
% - mean subtraction for the real + shuffles plot
% 20 September 2024
% - code clean-up

%% prep
%#ok<*AGROW>
latency = sLatenzy.latency;
peakTimes = sLatenzy.peakTimes;
peakVals = sLatenzy.peakVals;
realFrac = sLatenzy.realFrac;
fracLin = sLatenzy.fracLin;
realDiff = sLatenzy.realDiff;
realTime = sLatenzy.realTime;
meanRealDiff = sLatenzy.meanRealDiff;
randDiff = sLatenzy.randDiff;
randTime = sLatenzy.randTime;
meanRandDiff = sLatenzy.meanRandDiff;
peakZ = sLatenzy.peakZ;
latenzyIdx = sLatenzy.latencyIdx;

numIters = numel(peakTimes);

useColors = lines(numIters);
lineWidth = 1.5;
markerSize = 60; %50;

%% PLOT
figure;
figHandles = nan(1,4);
if makePlots==1
    %raster plot
    figHandles(1) = subplot(2,2,1);
    rasterPlot(spikeTimes,eventTimes,useMaxDur);
    hold on
    yLims = ylim;
    plot([latency latency], yLims, 'r--','LineWidth',lineWidth);
    set(gca,'box','off','TickDir','out');
    ylabel('Trial');
    xlabel('Time from event (s)');
    title('Spike raster plot');
end

%plot cumulative spike over time
figHandles(2) = subplot(2,2,2); hold on
for iter = 1:numIters
    plot(realTime{iter},fracLin{iter},'color',[0.5 0.5 0.5],'LineWidth',lineWidth);
    p(iter) =  plot(realTime{iter},realFrac{iter},'color',useColors(iter,:),'LineWidth',lineWidth);
    labels{iter} = sprintf('%d', iter);
end
set(gca,'box','off','TickDir','out');
xlim(useMaxDur);
xlabel('Time from event (s)');
ylabel('Fractional spike position');
title('Cumulative spikes');
lgd = legend(p, labels, 'Location', 'southeast','Box','off');
title(lgd, 'Iteration');

%plot offset from linear baseline
figHandles(3) = subplot(2,2,3); hold on
plot(useMaxDur,[0 0],'color',[0.5 0.5 0.5],'LineWidth',lineWidth,'LineStyle','--');
for iter = 1:numIters
    plot(realTime{iter},realDiff{iter},'color',useColors(iter,:),'LineWidth',lineWidth);
end
scatter(peakTimes(~latenzyIdx),peakVals(~latenzyIdx),markerSize,'xk','LineWidth',lineWidth);
scatter(peakTimes(latenzyIdx),peakVals(latenzyIdx),markerSize,'xr','LineWidth',lineWidth);
set(gca,'box','off','TickDir','out');
xlim(useMaxDur);
xlabel('Time from event (s)');
ylabel('Offset from lin. (Δfrac.)');
title('\color{red}x\color{black} = estimated onset');

%for iteration w/latency, plot real + shuffles
figHandles(4) = subplot(2,2,4); hold on
for thisShuffle=1:length(randDiff)
    plot(randTime{thisShuffle,latenzyIdx},(randDiff{thisShuffle,latenzyIdx}-meanRandDiff(thisShuffle,latenzyIdx)) ,'Color',[0.5 0.5 0.5 0.5],'LineWidth',lineWidth);
end
plot(realTime{latenzyIdx},(realDiff{latenzyIdx}-meanRealDiff(latenzyIdx)),'color',useColors(latenzyIdx,:),'LineWidth',lineWidth);
scatter(peakTimes(latenzyIdx),(peakVals(latenzyIdx)-meanRealDiff(latenzyIdx)),markerSize,'xr','LineWidth',lineWidth);
if latenzyIdx(1),xlim(useMaxDur);
else, xlim([useMaxDur(1)  peakTimes(find(latenzyIdx)-1)]);end
set(gca,'box','off','TickDir','out');
xlabel('Time from event (s)');
ylabel('Offset from lin. (Δfrac.)');
title(sprintf('Real data + shuffles, mean-subtracted (Z = %.1f)',peakZ(latenzyIdx)));

%add title
sgtitle(sprintf('Estimated response latency = %.4fs', latency), 'FontWeight', 'bold');
end