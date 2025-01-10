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
%   v0.9 - 6 January 2025
%   - created by Robin Haak

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
    sprintf('The second element of useMaxDur must be larger than the first element, you requested [%.3f %.3f]',...
    useMaxDur(1),useMaxDur(2)));

if ~exist('addArtifSpikes','var') || isempty(addArtifSpikes)
    addArtifSpikes = false;
end

%% compute relative spike times
spikesPerEvent = arrayfun(@(x) spikeTimes(spikeTimes > (x+useMaxDur(1)) & spikeTimes < (x+useMaxDur(2)))-x, ...
    eventTimes,'UniformOutput',false);

%concat and sort
relSpikeTimes = sort(vertcat(spikesPerEvent{:}));

%% if requested, add artificial spikes to cover full epoch
if addArtifSpikes && ~isempty(relSpikeTimes)
    relSpikeTimes = unique([useMaxDur(1); relSpikeTimes; useMaxDur(2)]);
end
end
