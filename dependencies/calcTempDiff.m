function [tempDiff,relSpikeTimes,spikeFrac,fracLinear] = calcTempDiff(spikeTimes,eventTimes,useMaxDur)
% compute temporal offset vector, syntax:
%   [tempDiff,relSpikeTimes,spikeFracs,fracLinear] = calcTempDiff(spikeTimes,eventTimes,useMaxDur)
%
% history:
%   v0.9 - 6 January 2025
%   - created by Robin Haak

%% prep
tempDiff = [];
spikeFrac = [];
fracLinear = [];

%get spikes relative to events (and add two artificial spikes)
relSpikeTimes = getRelSpikeTimes(spikeTimes,eventTimes,useMaxDur,true);
if isempty(relSpikeTimes)
    return
end

relSpikeTimes = getDistinctSpikes(relSpikeTimes);

%% get temporal offset vector
%fractional spike positions
numSpikes = numel(relSpikeTimes);
spikeFrac = linspace(1/numSpikes,1,numSpikes)';

%linear fraction
fracLinear = (relSpikeTimes-relSpikeTimes(1))./(relSpikeTimes(end)-relSpikeTimes(1));

%compute difference
tempDiff = spikeFrac-fracLinear;
% tempDiff = tempDiff-mean(tempDiff); % mean is subtracted in main latenzy() function

end