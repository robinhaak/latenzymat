function [pseudoSpikeTimes, pseudoEventTimes] = getPseudoTimes(spikeTimes,eventTimes,useMaxDur,discardEdges)
% perform data-stitching, syntax:
%   [pseudoSpikeTimes, pseudoEventTimes] = getPseudoTimes(spikeTimes,eventTimes,useMaxDur,discardEdges)
%
% history:
%   v0.9 - 6 January 2025
%   - created by Robin Haak

%ensure correct orientation
spikeTimes = sort(spikeTimes(:));
eventTimes = sort(eventTimes(:));

%check inputs
if ~exist('useMaxDur', 'var') || isempty(useMaxDur)
    eventTimes = sort(eventTimes);
    useMaxDur = min(diff(eventTimes));
end

if isscalar(useMaxDur), useMaxDur = sort([0 useMaxDur]); end
assert(useMaxDur(2) > useMaxDur(1), [mfilename ':WrongMaxDurInput'], ...
    sprintf('The second element of useMaxDur must be larger than the first element, you requested [%.3f %.3f]',...
    useMaxDur(1), useMaxDur(2)));

if ~exist('discardEdges', 'var') || isempty(discardEdges)
    discardEdges = false;
end

%initialize variables
sampleNum = numel(spikeTimes);
eventNum = numel(eventTimes);
pseudoSpikeT = cell(1, eventNum);
pseudoEventTimes = nan(eventNum, 1);
duration = useMaxDur(2)-useMaxDur(1);
pseudoEventT = 0;
lastUsedSample = 0;
firstSample = [];

%loop over each event
for thisEvent = 1:eventNum
    eventT = eventTimes(thisEvent)+useMaxDur(1);
    startSample = find(spikeTimes >= eventT, 1);
    endSample = find(spikeTimes > (eventT+duration), 1)-1;

    %handle cases where no pre-event spikes are found
    if isempty(startSample), startSample = sampleNum+1; end
    if isempty(endSample) || startSample > endSample
        startSample = [];
        endSample = sampleNum;
    end

    %handle edge case for first and last events without adding spikes
    theseSamples = startSample:endSample;
    remSamples = (theseSamples <= 0) | (theseSamples > sampleNum);
    useSamples = theseSamples(~remSamples);

    if ~isempty(useSamples)
        if thisEvent == 1 && ~discardEdges
            useSamples = 1:useSamples(end);
        elseif thisEvent == eventNum && ~discardEdges
            useSamples = useSamples(1):sampleNum;
        end
    end

    %handle spike times
    addT = spikeTimes(useSamples);
    samplesOverlap = (useSamples <= lastUsedSample);

    %manage overlaps between consecutive windows
    if thisEvent == 1
        pseudoEventT = 0;
    elseif thisEvent > 1 && duration > (eventT-eventTimes(thisEvent-1))
        useSamples = useSamples(~samplesOverlap);
        addT = spikeTimes(useSamples);
        pseudoEventT = pseudoEventT+eventT-eventTimes(thisEvent-1);
    else
        pseudoEventT = pseudoEventT+duration;
    end

    %make local pseudo times
    if isempty(useSamples)
        localPseudoT = [];
    else
        lastUsedSample = useSamples(end);
        localPseudoT = addT-eventT+pseudoEventT;
    end

    %set first sample for edge handling
    if isempty(firstSample) && ~isempty(useSamples)
        firstSample = useSamples(1);
        pseudoT0 = pseudoEventT;
    end

    pseudoSpikeT{thisEvent} = localPseudoT;
    pseudoEventTimes(thisEvent) = pseudoEventT;
end

%add beginning
if ~discardEdges && ~isempty(firstSample) && firstSample > 1
    stepBegin = spikeTimes(firstSample)-spikeTimes(firstSample-1);
    sampAddBeginning = 1:(firstSample-1);

    %only add beginning spikes if they exist
    if ~isempty(sampAddBeginning)
        pseudoSpikeT = cat(2,{spikeTimes(sampAddBeginning)-spikeTimes(sampAddBeginning(1))+...
            pseudoT0-stepBegin-range(spikeTimes(sampAddBeginning))},pseudoSpikeT);
    end
end

%add end
if ~discardEdges
    %find remaining spikes after the last event window
    lastUsedSample = find(spikeTimes > (eventTimes(end)+duration),1);
    if isempty(lastUsedSample)
        %if no spikes are beyond the last event window, take the remaining spikes
        sampAddEnd = (find(spikeTimes >= eventTimes(end),1)):sampleNum;
    else
        %otherwise, consider spikes within the range
        sampAddEnd = lastUsedSample:sampleNum;
    end

    %only add end spikes if they exist
    if ~isempty(sampAddEnd)
        pseudoSpikeT = cat(2,pseudoSpikeT,{spikeTimes(sampAddEnd)-eventTimes(end)+pseudoEventT});
    end
end

%recombine into a single vector
pseudoSpikeTimes = cell2mat(pseudoSpikeT(:));

%adjust event times
pseudoEventTimes = pseudoEventTimes+abs(useMaxDur(1));
end