function [tempDiff,relSpikeTimes,spikeFracs,fracLinear] = calcTempDiff2(spikeTimes,eventTimes1,eventTimes2,useMaxDur,useFastIn)
% compute temporal offset vector, syntax:
%   [tempDiff,relSpikeTimes,spikeFracs,fracLinear] = calcTempDiff(spikeTimes,eventTimes,useMaxDur)
%
% history:
%   v0.9 - 6 January 2025
%   - created by Robin Haak

%% prep
tempDiff = [];
spikeFracs = [];
fracLinear = [];

%get spikes relative to events (and add two artificial spikes)
relSpikeTimes = getRelSpikeTimes(spikeTimes,eventTimes,useMaxDur,true);
if isempty(relSpikeTimes)
    return
end

%introduce minimal jitter to repeating spike times (if any)
relSpikeTimes = sort(relSpikeTimes);
uniqueOffset = max(eps(relSpikeTimes));
idxRepeat = [false;diff(relSpikeTimes)<uniqueOffset];
while any(idxRepeat)
    notUnique = relSpikeTimes(idxRepeat);
    addJitter = cat(1,1+9*rand([numel(notUnique),1]),-1-9*rand([numel(notUnique),1]));
    addJitter = uniqueOffset*addJitter(randperm(numel(addJitter),numel(notUnique)));
    relSpikeTimes(idxRepeat) = relSpikeTimes(idxRepeat)+addJitter;
    relSpikeTimes = sort(relSpikeTimes);
    idxRepeat = [false;diff(relSpikeTimes)<uniqueOffset];
end

%% get temporal offset vector
%fractional spike positions
numSpikes = numel(relSpikeTimes);
spikeFracs = linspace(1/numSpikes,1,numSpikes)';

%linear fractions
fracLinear = (relSpikeTimes-relSpikeTimes(1))./(relSpikeTimes(end)-relSpikeTimes(1));

%compute difference
tempDiff = spikeFracs-fracLinear;
% tempDiff = tempDiff-mean(tempDiff); % mean is subtracted in main latenzy() function

end