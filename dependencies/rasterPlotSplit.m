function rasterPlotSplit(spikeTimes, eventTimes, useMaxDur, trialType, plotColor, plotMaxSpikes, doSort)
% make a raster plot, syntax:
%   rasterPlot(spikeTimes,eventTimes,useMaxDur,trialType,color,plotMaxSpikes,doSort)
%   input:
%   - spikeTimes [S x 1]: spike times (s)
%   - eventTimes [T x 1]: event (start) times (s)
%   - useMaxDur: scalar or [pre post], time to include after/around event times (s)
%   - trialType: label for trialtype (e.g., orientation), same size as eventTimes
%   - plotColor: colors for plotting [3 (RGB) x N], where N is unique trialtypes
%   - plotMaxSpikes: max number of spikes to plot (default: inf)
%   - doSort: boolean flag, sort spikes by trialType (default: true)
%
% history:
% 6 August 2024
%   - created by Robin Haak
% 19 September 2024
%   - code clean-up to make it more user friendly
%   - added new default colormap

%% prep 
%ensure correct orientation
spikeTimes = spikeTimes(:);
eventTimes = eventTimes(:);

%get inputs
if ~exist('useMaxDur', 'var') || isempty(useMaxDur)
    eventTimes = sort(eventTimes);
    useMaxDur = min(diff(eventTimes));
end
if isscalar(useMaxDur), useMaxDur = [0 useMaxDur]; end
assert(useMaxDur(1) <= 0, [mfilename ':WrongMaxDurInput'], ...
    sprintf('UseMaxDur(1) must be a negative scalar, you requested %.3f', useMaxDur(1)));
assert(useMaxDur(2) > 0, [mfilename ':WrongMaxDurInput'], ...
    sprintf('UseMaxDur(2) must be a positive scalar, you requested %.3f', useMaxDur(2)));

if ~exist('trialType', 'var') || isempty(trialType)
    trialType = ones(size(eventTimes));
end

%sort trial types and get unique values
[trialType, uniqueType, ~, ~, ~] = val2idx(trialType);
numTrialType = numel(uniqueType);

%ensure plotColor has the correct number of colors
if ~exist('plotColor', 'var') || isempty(plotColor) || size(plotColor, 1) ~= numTrialType
    colormapOut = generateCustomColorblindColormap(numTrialType-1);
    plotColor = [0 0 0; colormapOut];
end

%ensure plotMaxSpikes is defined
if ~exist('plotMaxSpikes', 'var') || isempty(plotMaxSpikes)
    plotMaxSpikes = inf;
end

%ensure doSort is defined
if ~exist('doSort', 'var') || isempty(doSort)
    doSort = true;
end

%subsample spikes if they exceed the maximum allowed for plotting
if numel(spikeTimes) > plotMaxSpikes
    spikeTimes = spikeTimes(sort(randperm(numel(spikeTimes), plotMaxSpikes)));
end

%% make raster plot
cla
hold all
offset = 0;
if doSort
    %sort by trial type
    for thisTrialType = 1:numel(uniqueType)
        %get event times for current trial type
        thisTrialStarts = eventTimes(trialType == thisTrialType);

        %get spike times in the specified window around events
        [~, spikesPerEvent] = getRelSpikeTimes(spikeTimes, thisTrialStarts, useMaxDur);

        %plot spikes for each trial in the current trial type
        for thisTrial = 1:numel(thisTrialStarts)
            vecTimes = spikesPerEvent{thisTrial};
            line([vecTimes(:)'; vecTimes(:)'], [(thisTrial + offset) * ones(1, numel(vecTimes)) - 0.5; (thisTrial + offset) * ones(1, numel(vecTimes)) + 0.5], ...
                'Color', plotColor(thisTrialType, :), 'LineWidth', 1.5);
        end
        offset = offset + numel(thisTrialStarts);  %update offset for next trial type
    end

else
    %don't sort by trial type
    [~, spikesPerEvent] = getRelSpikeTimes(spikeTimes, eventTimes, useMaxDur);

    %plot spikes for each trial
    for thisTrial = 1:numel(eventTimes)
        vecTimes = spikesPerEvent{thisTrial};
        line([vecTimes(:)'; vecTimes(:)'], [thisTrial * ones(1, numel(vecTimes)) - 0.5; thisTrial * ones(1, numel(vecTimes)) + 0.5], ...
            'Color', plotColor(trialType(thisTrial), :), 'LineWidth', 1.5);
    end
end
hold off

%set figure props
ylim([0.5 numel(eventTimes) + 0.5]);
xlim(useMaxDur);
xlabel('Time from event (s)');
ylabel('Trial');
end
