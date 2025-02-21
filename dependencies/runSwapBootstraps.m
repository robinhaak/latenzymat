function [peaksRandD,resampD,resampT] = runSwapBootstraps(spikesPerEvent1,spikesPerEvent2,useMaxDur,resampNum,useParPool,useFastInterp)
%run bootstraps by swapping trials, syntax:
%   [peaksRandD,resampD,resampT] = runSwapBootstraps(spikesPerEvent1,spikesPerEvent2,useMaxDur,resampNum,useParPool,useFastInterp)
%
% history:
%   v0.9 - 19 February 2025
%   - created by Robin Haak

%% run bootstraps
%swap trials randomly in each resampling
 %#ok<*PFBNS>
resampT = cell(1,resampNum);
resampD = cell(1,resampNum);
peaksRandD = nan(1,resampNum);
trialsAgg = cat(1,spikesPerEvent1,spikesPerEvent2);
idxSpikes0 = cellfun(@numel,trialsAgg)==0;
numEv1 = numel(spikesPerEvent1);
numEv2 = numel(spikesPerEvent2);
numEvTot = numEv1+numEv2;

if useParPool
    parfor resamp=1:resampNum
        %% get random subsample
        useRand1 = randi(numEvTot,[1,numEv1]);
        useRand2 = randi(numEvTot,[1,numEv2]);

        spikesPerTrial1_Rand = trialsAgg(useRand1);
        spikesPerTrial2_Rand = trialsAgg(useRand2);
        if all(idxSpikes0(useRand1)) && all(idxSpikes0(useRand2))
            continue;
        end

        %get difference
        [randD,randT] = calcTempDiff2(spikesPerTrial1_Rand,spikesPerTrial2_Rand,useMaxDur,useFastInterp);

        %get largest deviation
        maxVal = max(randD);
        minVal = min(randD);
        if abs(minVal) >= abs(maxVal)
            maxRandD = minVal;
        else
            maxRandD = maxVal;
        end

        resampD{resamp} = randD;
        resampT{resamp} = randT;
        if ~isnan(maxRandD)
            peaksRandD(resamp) = maxRandD;
        end
    end
else
    for resamp=1:resampNum
        %% get random subsample
        useRand1 = randi(numEvTot,[1,numEv1]);
        useRand2 = randi(numEvTot,[1,numEv2]);

        spikesPerTrial1_Rand = trialsAgg(useRand1);
        spikesPerTrial2_Rand = trialsAgg(useRand2);
        if all(idxSpikes0(useRand1)) && all(idxSpikes0(useRand2))
            continue;
        end

        %get difference
        [randD,randT] = calcTempDiff2(spikesPerTrial1_Rand,spikesPerTrial2_Rand,useMaxDur,useFastInterp);

        %get largest deviation
        maxVal = max(randD);
        minVal = min(randD);
        if abs(minVal) >= abs(maxVal)
            maxRandD = minVal;
        else
            maxRandD = maxVal;
        end

        resampD{resamp} = randD;
        resampT{resamp} = randT;
        if ~isnan(maxRandD)
            peaksRandD(resamp) = maxRandD;
        end
    end
end