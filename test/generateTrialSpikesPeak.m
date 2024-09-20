function [spikeTimes,eventTimes,onsetTimes,minLatency,useMaxDur,sponRates,peakRates,peakDurs,sustRates,jitterSize] = ...
    generateTrialSpikesPeak(numUnits,numTrials,minLatency,doJitter,jitterSize)
% generate trialwise spiking data with exponential ISI times
% Robin Haak, 2024
%
% NB, trial duration is hard-coded to 2s, with the onset latency varying
% from 0.1-1s (unless otherwise specified)

if ~exist('numUnits','var') || isempty(numUnits)
    numUnits = 100;
end

if ~exist('numTrials','var') || isempty(numTrials)
    numTrials = 20;
end

if ~exist('minLatency','var') || isempty(minLatency)
    minLatency = 0.1; %s
end

if ~exist('doJitter','var') || isempty(doJitter)
    doJitter = false;
end
if doJitter
    if ~exist('jitterSize','var') || isempty(jitterSize)
        jitterSize = 1e-2; %s
    end
else
    jitterSize = 0;
end

%% params
useMaxDur = 2; %s, trial duration
eventTimes = (3*useMaxDur):(2*useMaxDur):(2*useMaxDur*(numTrials+1));
stopTime = max(eventTimes)+4*useMaxDur;

%spontaneous rate
minSponRate = 1e-1;
maxSponRate = 1e2;
muSponRate = log((minSponRate+maxSponRate)/2);
sigSponRate =  (log(maxSponRate)-log(minSponRate))/4;
sponRates = lognrnd(muSponRate,sigSponRate,[1,numUnits*10]);
sponRates(sponRates<minSponRate | sponRates>maxSponRate) = [];
sponRates = randsample(sponRates,numUnits);

%evoked rate
minPeakRate = 1;
maxPeakRate = 1e2;
muPeakRate = log((minPeakRate+maxPeakRate)/2);
sigPeakRate = (log(maxPeakRate)-log(minPeakRate))/4;
peakRates = lognrnd(muPeakRate,sigPeakRate,[1,numUnits*10]);
peakRates(peakRates<minPeakRate | peakRates>maxPeakRate) = [];
peakRates = randsample(peakRates,numUnits);

%half of the units will get a sustained rate
sustRates = peakRates./(rand(1,numUnits)+2);
idxNoSust = randperm(numUnits,floor(numUnits/3));
sustRates(idxNoSust) = 0;

%response timing
onsetTimes = rand(1,numUnits) + minLatency;
minPeakDur = 5e-3;
maxPeakDur = 5e-1;
muPeakDur = log((minPeakDur+maxPeakDur)/2);
sigPeakDur = (log(maxPeakDur)-log(minPeakDur))/4;
peakDurs = lognrnd(muPeakDur,sigPeakDur,[1,numUnits*10]);
peakDurs(peakDurs<minPeakDur | peakDurs>maxPeakDur) = [];
peakDurs = randsample(peakDurs,numUnits);

%% generate spiking data
spikeTimes = cell(1,numUnits);
trialNum = numel(eventTimes);
for thisUnit = 1:numUnits
    thisSponR = sponRates(thisUnit);
    thisPeakR = peakRates(thisUnit);
    thisSustR = sustRates(thisUnit);
    thisOnset = onsetTimes(thisUnit);
    thisPeakD = peakDurs(thisUnit);

    if jitterSize > 0
        trialJitter = jitterSize*(2*rand(1,trialNum)-1);
    else
        trialJitter = 0;
    end

    epochRates = [0 thisPeakR thisSustR];
    epochDurs = [thisOnset thisPeakD];
    epochDurs(3) = useMaxDur-sum(epochDurs);

    sponSpikes = cumsum(exprnd(1/thisSponR,[1 round(thisSponR*stopTime*5)]));
    sponSpikes(sponSpikes>stopTime)=[];

    trialEpochSpikes = cell(numel(eventTimes),numel(epochDurs));
    for thisTrial=1:trialNum
        epochStart = eventTimes(thisTrial);

        if jitterSize > 0
            epochDurs = [thisOnset+trialJitter(thisTrial) thisPeakD];
            epochDurs(3) = useMaxDur-sum(epochDurs);
        end
        for thisEpoch=1:numel(epochDurs)
            thisDur = epochDurs(thisEpoch);
            thisEpochRate = epochRates(thisEpoch);
            thisEpochSpikes = cumsum(exprnd(1/thisEpochRate,[1 round(thisEpochRate*thisDur*5)]));
            thisEpochSpikes(thisEpochSpikes > thisDur) = [];
            trialEpochSpikes{thisTrial,thisEpoch} = epochStart+thisEpochSpikes;
            epochStart = epochStart + thisDur;
        end

    end
    spikeTimes{thisUnit} = sort(cat(1,sponSpikes(:),cell2vec(trialEpochSpikes)));
end

end