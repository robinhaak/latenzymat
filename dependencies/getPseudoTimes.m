function [pseudoSpikeTimes, pseudoEventTimes] = getPseudoTimes(spikeTimes, eventTimes, useMaxDur, discardEdges)
% perform data-stitching, syntax:
%   [pseudoSpikeTimes, pseudoEventTimes] = getPseudoTimes(spikeTimes, eventTimes, useMaxDur, discardEdges)

% Ensure correct orientation
spikeTimes = sort(spikeTimes(:));
eventTimes = sort(eventTimes(:));

% Check inputs
if ~exist('useMaxDur', 'var') || isempty(useMaxDur)
    eventTimes = sort(eventTimes);
    useMaxDur = min(diff(eventTimes));
end
if isscalar(useMaxDur), useMaxDur = [0 useMaxDur]; end
assert(useMaxDur(1) <= 0, sprintf('UseMaxDur(1) must be a negative scalar, you requested %.2f', useMaxDur(1)));
assert(useMaxDur(2) > 0, sprintf('UseMaxDur(2) must be a positive scalar, you requested %.2f', useMaxDur(2)));
if ~exist('discardEdges', 'var') || isempty(discardEdges)
    discardEdges = false;
end

% Initialize variables
sampleNum = numel(spikeTimes);
eventNum = numel(eventTimes);
pseudoSpikeT = cell(1, eventNum);
pseudoEventTimes = nan(eventNum, 1);
fullDuration = sum(abs(useMaxDur));
pseudoEventT = 0;
lastUsedSample = 0;
firstSample = [];

% Loop over each event
for thisEvent = 1:eventNum
    eventT = eventTimes(thisEvent) + useMaxDur(1);
    startSample = find(spikeTimes >= eventT, 1);
    endSample = find(spikeTimes > (eventT + fullDuration), 1) - 1;
    
    % Handle cases where no pre-event spikes are found
    if isempty(startSample), startSample = sampleNum + 1; end
    if isempty(endSample) || startSample > endSample
        startSample = [];
        endSample = [];
    end

    % Handle the case where no spikes exist in the current window
    if isempty(endSample)
        endSample = numel(spikeTimes);
    end
    
    theseSamples = startSample:endSample;
    remSamples = (theseSamples <= 0) | (theseSamples > sampleNum);
    useSamples = theseSamples(~remSamples);

    % Handle edge case for first and last events without adding spikes
    if ~isempty(useSamples)
        if thisEvent == 1 && ~discardEdges
            useSamples = 1:useSamples(end);
        elseif thisEvent == eventNum && ~discardEdges
            useSamples = useSamples(1):sampleNum;
        end
    end
    
    % Handle spike times
    addT = spikeTimes(useSamples);
    samplesOverlap = (useSamples <= lastUsedSample);

    % Manage overlaps between consecutive windows
    if thisEvent == 1
        pseudoEventT = 0;
    elseif thisEvent > 1 && fullDuration > (eventT - eventTimes(thisEvent - 1))
        useSamples = useSamples(~samplesOverlap);
        addT = spikeTimes(useSamples);
        pseudoEventT = pseudoEventT + eventT - eventTimes(thisEvent - 1);
    else
        pseudoEventT = pseudoEventT + fullDuration;
    end

    %% Make local pseudo times
    if isempty(useSamples)
        localPseudoT = [];
    else
        lastUsedSample = useSamples(end);
        localPseudoT = addT - eventT + pseudoEventT;
    end

    % Set first sample for edge handling
    if isempty(firstSample) && ~isempty(useSamples)
        firstSample = useSamples(1);
        pseudoT0 = pseudoEventT;
    end

    pseudoSpikeT{thisEvent} = localPseudoT;
    pseudoEventTimes(thisEvent) = pseudoEventT;
end

%% Add beginning
if ~discardEdges && ~isempty(firstSample) && firstSample > 1
    stepBegin = spikeTimes(firstSample) - spikeTimes(firstSample - 1);
    sampAddBeginning = 1:(firstSample - 1);
    
    % Only add beginning spikes if they exist
    if ~isempty(sampAddBeginning)
        pseudoSpikeT = cat(2, {spikeTimes(sampAddBeginning) - spikeTimes(sampAddBeginning(1)) + pseudoT0 - stepBegin - range(spikeTimes(sampAddBeginning))}, pseudoSpikeT);
    end
end

%% Add end
numSpikes = numel(spikeTimes);
lastUsedSample = find(spikeTimes > (eventTimes(end) + fullDuration), 1);

% Ensure no artificial spikes are added in the post-event window
if ~discardEdges && ~isempty(lastUsedSample) && numSpikes > lastUsedSample
    sampAddEnd = lastUsedSample:numSpikes;
    
    % Only add end spikes if they exist
    if ~isempty(sampAddEnd)
        pseudoSpikeT = cat(2, pseudoSpikeT, {spikeTimes(sampAddEnd) - eventT + pseudoEventT + fullDuration});
    end
end

%% Output
% Recombine into a single vector
pseudoSpikeTimes = cell2vec(pseudoSpikeT);

% Adjust event times
pseudoEventTimes = pseudoEventTimes + abs(useMaxDur(1));
end

% function [pseudoSpikeTimes,pseudoEventTimes] = getPseudoTimes(spikeTimes,eventTimes,useMaxDur,discardEdges)
% % perform data-stitching, syntax:
% %   [pseudoSpikeTimes,pseudoEventTimes] = getPseudoTimes(spikeTimes,eventTimes,useMaxDur,discardEdges)
% %
% % history:
% % 5 Feb 2024
% %   - created by Robin Haak based almost entirely on getPseudoSpikeVectors() v3.5.2 (a zetatest() subfunction by Jorrit Montijn),
% %       main difference is handling of pre-event time)
% % 7 Feb 2024
% %   - minor changes to the code
% % 10 April 2024
% %   - useMaxDur instead of separate pre- and post-event time variables
% 
% %% prep
% %ensure correct orientation
% spikeTimes = sort(spikeTimes(:));
% eventTimes = sort(eventTimes(:));
% 
% %check inputs
% if ~exist('useMaxDur','var') || isempty(useMaxDur)
%     eventTimes = sort(eventTimes);
%     useMaxDur = min(diff(eventTimes));
% end
% if isscalar(useMaxDur), useMaxDur = [0 useMaxDur]; end
% assert(useMaxDur(1)<=0,[mfilename ':WrongMaxDurInput'],...
%     sprintf('UseMaxDur(1) must be a negative scalar, you requested %.2f',useMaxDur(1)));
% assert(useMaxDur(2)>0,[mfilename ':WrongMaxDurInput'],...
%     sprintf('UseMaxDur(2) must be a positive scalar, you requested %.2f',useMaxDur(2)));
% 
% if ~exist('discardEdges','var') || isempty(discardEdges)
%     discardEdges = false;
% end
% 
% %% run
% sampleNum = numel(spikeTimes);
% eventNum = numel(eventTimes);
% pseudoSpikeT = cell(1,eventNum);
% pseudoEventTimes = nan(eventNum,1);
% fullDuration = sum(abs(useMaxDur));
% pseudoEventT = 0;
% lastUsedSample = 0;
% firstSample = [];
% for thisEvent = 1:eventNum
%     eventT = eventTimes(thisEvent)+useMaxDur(1);
%     startSample = (find(spikeTimes >= eventT,1));
%     endSample = (find(spikeTimes > (eventT+fullDuration),1)-1);
%     if startSample > endSample
%         endSample = [];
%         startSample = [];
%     end
%     if isempty(endSample)
%         endSample = numel(spikeTimes);
%     end
%     theseSamples = startSample:endSample;
%     remSamples = (theseSamples<=0) | (theseSamples>sampleNum);
%     useSamples = theseSamples(~remSamples);
% 
%     %check if beginning or end
%     if ~isempty(useSamples)
%         if thisEvent==1 && ~discardEdges
%             useSamples = 1:useSamples(end);
%         elseif thisEvent==eventNum && ~discardEdges
%             useSamples = useSamples(1):sampleNum;
%         end
%     end
%     addT = spikeTimes(useSamples);
%     samplesOverlap = (useSamples<=lastUsedSample);
% 
%     %get event t
%     if thisEvent == 1
%         pseudoEventT = 0;
%     else
%         if thisEvent > 1 && fullDuration > (eventT-eventTimes(thisEvent-1))
%             useSamples = useSamples(~samplesOverlap);
%             addT = spikeTimes(useSamples);
%             pseudoEventT = pseudoEventT+eventT-eventTimes(thisEvent-1);
%         else
%             pseudoEventT = pseudoEventT+fullDuration;
%         end
%     end
% 
%     %% MAKE LOCAL TO EVENT
%     if isempty(useSamples)
%         localPseudoT = [];
%     else
%         lastUsedSample = useSamples(end);
%         localPseudoT = addT-eventT+pseudoEventT;
%     end
%     if isempty(firstSample) && ~isempty(useSamples)
%         firstSample = useSamples(1);
%         pseudoT0 = pseudoEventT;
%     end
% 
%     pseudoSpikeT{thisEvent} = localPseudoT;
%     pseudoEventTimes(thisEvent) = pseudoEventT;
% 
% end
% 
% %% add beginning
% if ~discardEdges && ~isempty(firstSample) && firstSample > 1
%     stepBegin = spikeTimes(firstSample)-spikeTimes(firstSample-1);
%     sampAddBeginning = 1:(firstSample-1);
%     pseudoSpikeT = cat(2,{spikeTimes(sampAddBeginning)-...
%         spikeTimes(sampAddBeginning(1))+pseudoT0-stepBegin-range(spikeTimes(sampAddBeginning))},pseudoSpikeT);
% end
% 
% %% add end
% numSpikes = numel(spikeTimes);
% lastUsedSample = find(spikeTimes>(eventTimes(end)+fullDuration),1);
% if ~discardEdges && ~isempty(lastUsedSample) && numSpikes > lastUsedSample
%     sampAddEnd = lastUsedSample:numSpikes;
%     pseudoSpikeT = cat(2,pseudoSpikeT,{spikeTimes(sampAddEnd)-eventT+pseudoEventT+fullDuration});
% end
% 
% %% out
% %recombine into vector
% pseudoSpikeTimes = cell2vec(pseudoSpikeT);
% 
% %re-adjust event times
% pseudoEventTimes = pseudoEventTimes+abs(useMaxDur(1));
% end
