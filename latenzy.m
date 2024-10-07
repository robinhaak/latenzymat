function [respLatency,sLatenzy] = latenzy(spikeTimes,eventTimes,useMaxDur,resampNum,jitterSize,minPeakZ,doStitch,useParPool,useDirectQuant,allowNegative,makePlots)
% get stimulus-associated response latency, syntax:
%   [respLatency,sLatenzy] = latenzy(spikeTimes,eventTimes,useMaxDur,resampNum,jitterSize,minPeakZ,doStitch,useParPool,useDirectQuant,allowNegative,makePlots)
%   input:
%   - spikeTimes: [S x 1]: spike times (s)
%   - eventTimes: [T x 1]: event (start) times (s)
%   - useMaxDur: scalar or [N x 2], time to include after/around event times (s) (default: [0 min(diff(eventtimes))])
%   - resampNum: integer, number of resamples (default: 100)
%   - jitterSize: scalar, temporal jitter window relative to useMaxDur (s) (default: 2)
%   - minPeakZ: scalar, minimal z-score for peak significance (default: 2)
%   - doStitch: boolean flag, perform data stitching, highly recommended! (default: true)
%   - useParPool: boolean flag, use parallel pool for resamples (default: true, but only when parallel pool is already active!)
%   - useDirectQuant: boolean flag, use the empirical null-distribution rather than the Gumbel approximation (default: false)
%   - allowNegative: boolean flag, allow negative latencies (default: false)
%   - makePlots: integer, plotting switch (0=none, 1=raster+traces, 2=traces only, default: 0)
%
%   output:
%   - respLatency: response latency (s)
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
%       - pValsPeak: p-values obtained by comparing to shuffles
%       - peakZ: same as pValsPeak but z-scored
%       - latencyIdx: use to index arrays above
%       - handleFigs: figure handles
%
% history:
% 11 April 2024
% - created by Robin Haak
% 23 April 2024
% - cleaned up the code
% 8 August 2024
% - worked on plotting function
% - changed peakSwitch so that it is always 0
% 12 August 2024
% - worked on computation of significance
% 16 September 2024
% - mean subtraction to make computation of significance time-invariant (https://elifesciences.org/articles/71969, see 'a proof of time-invariance')
% - removed restrictRange as input
% - changed parallel pool behavior: parfor enabled by default when parpool is active
% 18 September 2024
% - added allowNegative flag, to allow negative latencies (e.g., for behavioral events)
% 20 September 2024
% - made sure that code runs when Parallel Computing Toolbox is not present

%% prep
%ensure correct orientation
spikeTimes = spikeTimes(:);
eventTimes = eventTimes(:);

%get useMaxDur
if ~exist('useMaxDur','var') || isempty(useMaxDur)
    eventTimes = sort(eventTimes);
    useMaxDur = min(diff(eventTimes));
end
if isscalar(useMaxDur),useMaxDur = [0 useMaxDur];end
assert(useMaxDur(1)<=0,[mfilename ':WrongMaxDurInput'],...
    sprintf('The first element of useMaxDur must be a negative scalar, you requested %.3f',useMaxDur(1)));
assert(useMaxDur(2)>0,[mfilename ':WrongMaxDurInput'],...
    sprintf('The second element of useMaxDur must be a positive scalar, you requested %.3f',useMaxDur(2)));

%get resampNum
if ~exist('resampNum','var') || isempty(resampNum)
    resampNum = 100;
end

%get jitterSize
if ~exist('jitterSize','var') || isempty(jitterSize)
    jitterSize = 2;
end
assert(jitterSize>0,[mfilename ':WrongJitterInput'], ...
    sprintf('jitterSize must be a positive scalar, you requested %.3f',jitterSize));

%get peakZ
if ~exist('minPeakZ','var') || isempty(minPeakZ)
    minPeakZ = 2;
end

%get doStitch
if ~exist('doStitch','var') || isempty(doStitch)
    doStitch = true;
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
if ~exist('allowNegative','var') || isempty(allowNegative)
    allowNegative = false;
end

%get makePlots
if ~exist('makePlots','var') || isempty(makePlots)
    makePlots = 0;
end

%% MAIN
%pre-allocate
%#ok<*AGROW>
respLatency = nan;
peakTimesAgg = [];
peakValsAgg = [];
realFracAgg = {};
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
sLatenzy = struct();

%check if negative latencies are allowed
minLatency = 0;
if allowNegative
    minLatency = useMaxDur(1);
end

%run
while doContinue
    thisIter = thisIter+1;

    %perform data-stitching
    if doStitch
        discardEdges = true;
        [pseudoSpikeTimes,pseudoEventTimes] = getPseudoTimes(spikeTimes,eventTimes,thisMaxDur,discardEdges);
    else
        pseudoSpikeTimes = spikeTimes;
        pseudoEventTimes = eventTimes;
    end

    %get largest deviation
    [realDiff,realTime,spikeFracs,fracLinear] = calcTempOffset(pseudoSpikeTimes,pseudoEventTimes,thisMaxDur);
    [realPeakV,realPeakT] = findPeak(realDiff,realTime);
    realPeakSub = realPeakV-mean(realDiff);

    %run bootstraps
    [peaksRand,randDiff,randTime] = jitteredBootstraps(pseudoSpikeTimes,pseudoEventTimes,...
        thisMaxDur,resampNum,jitterSize,0,useParPool);
    meanRandDiff = cellfun(@(x)mean(x),randDiff);
    peaksRandSub = peaksRand-meanRandDiff;

    %compute significance
    [pValPeak,peakZ] = getPvalZ(abs(realPeakSub),peaksRandSub(~isnan(peaksRandSub)),useDirectQuant);

    %store
    if ~isnan(realPeakT)
        peakValsAgg(1,thisIter) = realPeakV;
        peakTimesAgg(1,thisIter) = realPeakT;
        realFracAgg{1,thisIter} = spikeFracs;
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
    if realPeakT > minLatency && peakZ>minPeakZ && ~isinf(peakZ)
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
else
    return
end

%build output
sLatenzy = struct;
sLatenzy.latency = respLatency;
sLatenzy.peakTimes = peakTimesAgg;
sLatenzy.peakVals = peakValsAgg;
sLatenzy.realFrac = realFracAgg;
sLatenzy.fracLin = fracLinAgg;
sLatenzy.realDiff = realDiffAgg;
sLatenzy.realTime = realTimeAgg;
sLatenzy.meanRealDiff = meanRealDiffAgg;
sLatenzy.randDiff = randDiffAgg;
sLatenzy.randTime = randTimeAgg;
sLatenzy.meanRandDiff = meanRandDiffAgg;
sLatenzy.pValsPeak = pValPeakAgg;
sLatenzy.peakZ = peakZAgg;
sLatenzy.latencyIdx = peakTimesAgg==respLatency;

%plot, optional
if makePlots>0, sLatenzy.figHandles = makeLatenzyFigs(sLatenzy,spikeTimes,eventTimes,useMaxDur,makePlots); end

end