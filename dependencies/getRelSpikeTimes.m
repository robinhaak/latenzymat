function [relSpikeTimes,spikesPerEvent] = getRelSpikeTimes(spikeTimes,eventTimes,useMaxDur,addArtifSpikes)
% create a vector of spike times relative to event times, syntax:
%   [relSpikeTimes,spikesPerEvent] = getRelSpikeTimes(spikeTimes,eventTimes,useMaxDur,addArtifSpks)
%   input:
%   - spikeTimes [S x 1]: spike times (s)
%   - eventTimes [T x 1]: event (start) times (s)
%   - useMaxDur: scalar or [pre post], time to include after/around event times (s)
%   - addArtifSpikes: boolean, add artificial spikes at beginning and end of epoch (default: true)
%
%   output:
%   - relSpikeTimes: spike times relative to events (s), sorted
%   - spikesPerEvent: relative spike times per event (s), sorted
%
% history:
% 1 Aug 2023
%   - created by Robin Haak
% 28 Aug 2023
%   - added spikesPerEvent as (optional) output
% 5 Feb 2024
%   - minor changes to code
% 9 April 2024
%   - added option to add two artificial spikes, at beginning and end
%   - useMaxDur instead of separate pre- and post-event time variables
% 14 November 2024
%   - changed useMaxDur behavior

%% prep
%ensure correct orientation
spikeTimes = spikeTimes(:);
eventTimes = eventTimes(:);

%check inputs
if ~exist('useMaxDur','var') || isempty(useMaxDur)
    eventTimes = sort(eventTimes);
    useMaxDur = min(diff(eventTimes));
end

if isscalar(useMaxDur), useMaxDur = sort([0 useMaxDur]); end
assert(useMaxDur(2)>useMaxDur(1),[mfilename ':WrongMaxDurInput'],...
    sprintf('The second element of useMaxDur must be larger than the first element, you requested [%.3f %.3f]',useMaxDur(1),useMaxDur(2)));

if ~exist('addArtifSpikes','var') || isempty(addArtifSpikes)
    addArtifSpikes = false;
end

%% compute relative spike times
numEvents = length(eventTimes);
spikesPerEvent = cell(numEvents,1);
for event = 1:numEvents
    startTime = eventTimes(event);
    stopTime = startTime+useMaxDur(end);
    spikesPerEvent{event} = ...
        spikeTimes(spikeTimes>(startTime+useMaxDur(1)) & spikeTimes<stopTime)-startTime;
end

%sort ascending
spikesPerEvent = cellfun(@(x)sort(x),spikesPerEvent,'UniformOutput',false);
relSpikeTimes = vertcat(spikesPerEvent{:});
relSpikeTimes = sort(relSpikeTimes);

%% if requested, add artificial spikes to cover full epoch
if addArtifSpikes && ~isempty(relSpikeTimes)
    if relSpikeTimes(1) > useMaxDur(1)
        relSpikeTimes = [useMaxDur(1); relSpikeTimes];
    end
    if relSpikeTimes(end) < useMaxDur(end)
        relSpikeTimes = [relSpikeTimes; useMaxDur(2)];
    end
end
end