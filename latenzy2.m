function [respLatency,sLatenzy2] = latenzy2(spikeTimes1,eventTimes1,spikeTimes2,eventTimes2,useMaxDur,resampNum,minPeakZ,useParPool,useDirectQuant,restrictNeg,makePlots)
% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx, syntax:
% [respLatency,sLatenzy] = latenzy2(spikeTimes1,eventTimes1,spikeTimes2,eventTimes2,useMaxDur,resampNum,minPeakZ,useParPool,useDirectQuant,restrictNeg,makePlots)  
%   inputs:
%   - spikeTimes: [S x 1]: spike times (s)
%   - eventTimes: [T x 1]: event (start) times (s)
%   - useMaxDur: scalar or [N x 2], time to include after/around event times (s) (default: [0 min(diff(eventtimes))])
%   - resampNum: integer, number of resamples (default: 100)
%   - jitterSize: scalar, temporal jitter window relative to useMaxDur (s) (default: 2)
%   - minPeakZ: scalar, minimal z-score for significance (default: 1.96)
%   - doStitch: boolean flag, perform data stitching, highly recommended! (default: true)
%   - useParPool: boolean flag, use parallel pool for resamples (default: true, but only when parallel pool is already active!)
%   - useDirectQuant: boolean flag, use the empirical null-distribution rather than the Gumbel approximation (default: false)
%   - restrictNeg: boolean flag, restrict negative latencies (default: true)
%   - makePlots: integer, plotting switch (0=none, 1=raster+traces, 2=traces only, default: 0)
%
%   outputs:
%   - respLatency: response latency (s) (NaN when no latency could be estimated)
%   - sLatenzy: structure with fields:
%       - latency: response latency (s)
%       - peakTimes: detected peak/through times, one per iter (s)
%       - peakVals: detected peak/through values, one per iter
%       - realFrac: used for plotting only
%       - fracLin: idem
%       - realDiff: idem
%       - realTime: idem
%       - meanRealDiff: idem
%       - randDiff: idem
%       - randTime: idem
%       - meanRandDiff: idem
%       - peakZ: significance z-scores
%       - pValsPeak: p-values corresponding to peakZ
%       - latenzyIdx: use to index arrays above
%       - handleFigs: figure handles
%
% history:
%   v0.9 - 6 January 2025
%   - created by Robin Haak

%% prep
%ensure correct orientation
spikeTimes1 = spikeTimes1(:);
if isempty(spikeTimes2)
    spikeTimes2 = spikeTimes1;
else
    spikeTimes2 = spikeTimes2(:);
end
eventTimes1 = eventTimes1(:);
eventTimes2 = eventTimes2(:);

%get useMaxDur
if ~exist('useMaxDur','var') || isempty(useMaxDur)
    eventTimes1 = sort(eventTimes1);
    eventTimes2 = sort(eventTimes2);
    useMaxDur = min([min(diff(eventTimes1)) min(diff(eventTimes2))]);
end
if isscalar(useMaxDur)
    useMaxDur = sort([0 useMaxDur]);
elseif numel(useMaxDur)~=2
    error([mfilename ':WrongMaxDurInput'],'useMaxDur must be a scalar or a two-element array');
end

%check useMaxDur
if useMaxDur(2)>0
    assert(useMaxDur(1)<=0,[mfilename ':WrongMaxDurInput'],...
        sprintf('When useMaxDur(2) > 0, useMaxDur(1) must be a negative scalar or 0, you requested [%.3f %.3f]',useMaxDur(1),useMaxDur(2)));
elseif useMaxDur(2)==0
    assert(useMaxDur(1)<0,[mfilename ':WrongMaxDurInput'],...
        sprintf('When useMaxDur(2) is 0, useMaxDur(1) must be a negative scalar, you requested [%.3f %.3f]',useMaxDur(1),useMaxDur(2)));
elseif useMaxDur(2)<0
        error([mfilename ':WrongMaxDurInput'],'useMaxDur(2) cannot be negative when useMaxDur(1) is negative!');
end

%get resampNum
if ~exist('resampNum','var') || isempty(resampNum)
    resampNum = 100; %zeta does 250?
end

%get peakZ
if ~exist('minPeakZ','var') || isempty(minPeakZ)
    minPeakZ = 1.96; %corresponds to alpha=0.05
end

%get useParPool
if ~exist('useParPool','var') || isempty(useParPool)
    try
        objPool = gcp('nocreate'); %get current parallel pool (no creation)
        if isempty(objPool) || ~isprop(objPool,'NumWorkers') || objPool.NumWorkers < 4
            useParPool = false; 
        else
            useParPool = true;
        end
    catch
        useParPool = false;
    end
end

%useDirectQuant
if ~exist('useDirectQuant','var') || isempty(useDirectQuant)
    useDirectQuant = false;
end

%allowNegative
if ~exist('restrictNegative','var') || isempty(restrictNeg)
    restrictNeg = true;
end

%get makePlots
if ~exist('makePlots','var') || isempty(makePlots)
    makePlots = 0;
end

%enable warning for late
giveLateWarn = true;

%% MAIN
%pre-allocate
%#ok<*AGROW>
respLatency = nan;
peakTimesAgg = [];
peakValsAgg = [];
realFracAgg = {};
diffUnSubAgg = {};
fracLinAgg = {};
realDiffAgg = {};
realTimeAgg = {};
meanRealDiffAgg = [];
randDiffAgg = {};
randTimeAgg = {};
meanRandDiffAgg = [];
pValPeakAgg = [];
peakZAgg = [];
keepPeaks = [];
thisMaxDur = useMaxDur;
doContinue = true;
thisIter = 0;
sLatenzy2 = struct;

%check if negative latencies are restricted
minLatency = useMaxDur(1);
if restrictNeg
    minLatency = 0;
end

% check if we can use the fast interpolation
useFastInterp = false;
try
    vecTest = lininterp1f([0;0.25;1],[0;0.5;1],[-1;0;0.1;0.5;1],nan);
    if isnan(vecTest(1)) && all(vecTest(2:5)==[0 0.2 2/3 1])
        useFastInterp = true;
    else
        error('lininterp1f gives incorrect output');
    end
catch
    try
        mex('lininterp1f.c');
        vecTest = lininterp1f([0;0.25;1],[0;0.5;1],[-1;0;0.1;0.5;1],nan);
        if isnan(vecTest(1)) && all(vecTest(2:5)==[0 0.2 2/3 1])
            useFastInterp = true;
        else
            error('lininterp1f gives incorrect output');
        end
    catch
        useParPool = false;
    end
end

%%
%run
while doContinue
    thisIter = thisIter+1;

    %get spikes per event
    [~,spikesPerEvent1] = getRelSpikeTimes(spikeTimes1,eventTimes1,thisMaxDur);
    [~,spikesPerEvent2] = getRelSpikeTimes(spikeTimes2,eventTimes2,thisMaxDur);

    %get temporal difference
    [realDiff,realTime,spikeFrac1,~,spikeFrac2,~,tempDiffUnSub,fracLinear] = ...
        calcTempDiff2(spikesPerEvent1,spikesPerEvent2,thisMaxDur,useFastInterp);
    if numel(realDiff) < 3
        return
    end

    %get largest deviation
    [maxDiff,maxIdx] = max(realDiff);
    [minDiff,minIdx] = min(realDiff);
    if abs(minDiff) >= abs(maxDiff)
        realMaxD = minDiff;
        realMaxIdx = minIdx;
    else
        realMaxD = maxDiff;
        realMaxIdx = maxIdx;
    end
    realPeakT = realTime(realMaxIdx);
    realPeakSub = realMaxD-mean(realDiff);

    %run bootstraps
    [peaksRand,randDiff,randTime] = runSwapBootstraps(spikesPerEvent1,spikesPerEvent2,...
        useMaxDur,resampNum,useParPool,useFastInterp);
    meanRandDiff = cellfun(@(x)mean(x),randDiff);
    peaksRandSub = peaksRand-meanRandDiff;

    %compute significance
    [pValPeak,peakZ] = computeZ(abs(realPeakSub),peaksRandSub(~isnan(peaksRandSub)),useDirectQuant);

    %store
    if ~isnan(realPeakT)
        peakValsAgg(1,thisIter) = realMaxD;
        peakTimesAgg(1,thisIter) = realPeakT;
        realFracAgg{1,thisIter} = spikeFrac1;
        realFracAgg{2,thisIter} = spikeFrac2;
        diffUnSubAgg{1,thisIter} = tempDiffUnSub;
        fracLinAgg{1,thisIter} = fracLinear;
        realDiffAgg{1,thisIter} = realDiff;
        realTimeAgg{1,thisIter} = realTime;
        meanRealDiffAgg(1,thisIter) = mean(realDiff);
        randDiffAgg(:,thisIter) = randDiff;
        randTimeAgg(:,thisIter) = randTime;
        meanRandDiffAgg(:,thisIter) = meanRandDiff;
        pValPeakAgg(1,thisIter) = pValPeak;
        peakZAgg(1,thisIter) = abs(peakZ);
    end

    %check whether to continue
    if realPeakT > minLatency && peakZ > minPeakZ && ~isinf(peakZ)
        keepPeaks(thisIter) = true;
        thisMaxDur(2) = realPeakT;
    else
        doContinue = false;
        keepPeaks(thisIter) = false;
    end
end

%get latency
thesePeakTimes = peakTimesAgg(logical(keepPeaks));
if ~isempty(thesePeakTimes)
    respLatency = thesePeakTimes(end);

    %warning
    if respLatency > (useMaxDur(1)+sum(abs(useMaxDur))/2) && giveLateWarn
        warning('Estimated latency is quite late in the window (>T/2), consider plotting and/or adjusting window');
    end
else
    return
end

%build output
sLatenzy2.latency = respLatency;
sLatenzy2.peakTimes = peakTimesAgg;
sLatenzy2.peakVals = peakValsAgg;
sLatenzy2.realFrac = realFracAgg;
sLatenzy2.diffUnSub = diffUnSubAgg;
sLatenzy2.fracLin = fracLinAgg;
sLatenzy2.realDiff = realDiffAgg;
sLatenzy2.realTime = realTimeAgg;
sLatenzy2.meanRealDiff = meanRealDiffAgg;
sLatenzy2.randDiff = randDiffAgg;
sLatenzy2.randTime = randTimeAgg;
sLatenzy2.meanRandDiff = meanRandDiffAgg;
sLatenzy2.pValsPeak = pValPeakAgg;
sLatenzy2.peakZ = peakZAgg;
sLatenzy2.latenzyIdx = peakTimesAgg==respLatency;

%plot, optional
if makePlots>0, sLatenzy2.figHandles = makeLatenzy2Figs(sLatenzy2,spikeTimes1,eventTimes1,spikeTimes2,eventTimes2,useMaxDur,makePlots); end

end