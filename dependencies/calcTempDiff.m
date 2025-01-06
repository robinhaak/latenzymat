function [tempDiff,relSpikeTimes,spikeFracs,fracLinear] = calcTempDiff(spikeTimes,eventTimes,useMaxDur)
% compute temporal offset vector, syntax:
%   [tempDiff,relSpikeTimes,spikeFracs,fracLinear] = calcTempDiff(spikeTimes,eventTimes,useMaxDur)
%
% history:
% 2 Aug 2023
%   - created by Robin Haak
% 7 Aug 2023
%   - minor changes to code
% 9 Aug 2023
%   - added jittering of repeating spike times
% 28 Aug 2023
%   - getRelSpikeTimes() is now called from within this function
% 5 Feb 2024
%   - updated jittering of repeating spiketimes, making it identical to getUniqueSpikes() v3.1.0 (a zetatest() subfunction by Jorrit Montijn)
% 10 April 2024
%   - useMaxDur instead of separate pre- and post-event time variables
%   - changed computation of offset vector so that it now gives the deviation in fractions instead of spike counts

%% prep
tempDiff = [];
spikeFracs = [];
fracLinear = [];

%get spikes relative to events
relSpikeTimes = getRelSpikeTimes(spikeTimes,eventTimes,useMaxDur,true);
if isempty(relSpikeTimes)
    return
end

%introduce minimal jitter to repeating spike times (if any)
relSpikeTimes = sort(relSpikeTimes);
uniqueOffset = max(eps(relSpikeTimes));
indRepeat = [false;diff(relSpikeTimes)<uniqueOffset];
while any(indRepeat)
    notUnique = relSpikeTimes(indRepeat);
    addJitter = cat(1,1+9*rand([numel(notUnique),1]),-1-9*rand([numel(notUnique),1]));
    addJitter = uniqueOffset*addJitter(randperm(numel(addJitter),numel(notUnique)));
    relSpikeTimes(indRepeat) = relSpikeTimes(indRepeat)+addJitter;
    relSpikeTimes = sort(relSpikeTimes);
    indRepeat = [false;diff(relSpikeTimes)<uniqueOffset];
end

%% get temporal offset vector
%fractional spike positions
spikeFracs = linspace(1/numel(relSpikeTimes),1,numel(relSpikeTimes))';

%linear fractions
fracLinear = (relSpikeTimes-relSpikeTimes(1))./(relSpikeTimes(end)-relSpikeTimes(1));

%compute difference
tempDiff = spikeFracs-fracLinear;
% tempDiff = tempDiff-mean(tempDiff); mean is subtracted in main latenzy() function
end