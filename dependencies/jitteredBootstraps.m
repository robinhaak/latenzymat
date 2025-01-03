function [peaksRandD,resampD,resampT] = jitteredBootstraps(spikeTimes,eventTimes,useMaxDur,resampNum,jitterSize,switchPosNeg,useParPool)
% run bootstraps by jittering event times, syntax:
%   [peaksRandD,randD,randT] = jitteredBootstraps(spikeTimes,eventTimes,preEventTime,postEventTime,resampNum,jitterSize,switchPosNeg,useParPool)
%
% history:
% 9 Feb 2024
%   - created by Robin Haak
% 15 Feb 2024
%   - cleanup of code
% 10 April 2024
%   - useMaxDur instead of separate pre- and post-event time variables
% 30 December 2024
% - removed separate function for peak detection, using min/max instead

%% prep
%ensure correct orientation
spikeTimes = sort(spikeTimes(:));
eventTimes = sort(eventTimes(:));

%% create random bootstraps
resampD = cell(1,resampNum);
resampT = cell(1,resampNum);
peaksRandD = nan(1,resampNum);
eventNum = numel(eventTimes);
fullDuration = diff(useMaxDur);
jitterPerTrial = nan(eventNum,resampNum);

%uniform jitters between jitterSize*[-tau,+tau]
for resamp=1:resampNum
    jitterPerTrial(:,resamp) = jitterSize*fullDuration*((rand(size(eventTimes))-0.5)*2);
end

%% run bootstraps
if useParPool
    parfor resamp=1:resampNum
        randEventT = eventTimes+jitterPerTrial(:,resamp);
        [randD,randT] = calcTempOffset(spikeTimes,randEventT,useMaxDur);
        % peakRandD = findPeak(randD,randT,[],switchPosNeg);
    
        maxVal = max(randD);
        minVal = min(randD);
        if abs(minVal) >= abs(maxVal)
            peakRandD = minVal;
        else
            peakRandD = maxVal;
        end

        resampD{resamp} = randD;
        resampT{resamp} = randT;
        if ~isnan(peakRandD)
            peaksRandD(resamp) = peakRandD;
        end
    end
else
    for resamp=1:resampNum
        randEventT = eventTimes+jitterPerTrial(:,resamp);
        [randD,randT] = calcTempOffset(spikeTimes,randEventT,useMaxDur);
        % peakRandD = findPeak(randD,randT,[],switchPosNeg);
 
        maxVal = max(randD);
        minVal = min(randD);
        if abs(minVal) >= abs(maxVal)
            peakRandD = minVal;
        else
            peakRandD = maxVal;
        end
        resampD{resamp} = randD;
        resampT{resamp} = randT;
        if ~isnan(peakRandD)
            peaksRandD(resamp) = peakRandD;
        end
    end
end

end