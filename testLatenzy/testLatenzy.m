% testLatenzy
%
% To run use runtests('testLatenzy')
%
% 2025, Robin Haak, Alexander Heimel

pth = fileparts(which('testLatenzy'));
load(fullfile(pth,'Topo2_20220126_AP.mat'));
idxIncl = ([sAP.sCluster.KilosortGood] | [sAP.sCluster.Contamination] < 0.1) & ...
    contains({sAP.sCluster.Area},'primary visual','IgnoreCase',true);
spikeTimesAgg = {sAP.sCluster(idxIncl).SpikeTimes};
sBlock = sAP.cellBlock{4}; %DG, 1s
eventTimes = sBlock.vecStimOnTime;
useParPool = false;
selClus = [16 19 20];
numClus = length(selClus);
res = NaN(1,3 * numClus); % 1 x n_tests * n_clus
saved = load(fullfile(pth,'testResults.mat'),'res');

%% Test 1, defaults
count = 1;
i = 1;
thisClus = selClus(i);
theseSpikeTimes = spikeTimesAgg{thisClus};

useMaxDur = 1;
rng(1,'twister');
res(count) = latenzy(theseSpikeTimes,eventTimes,useMaxDur,[],[],[],[],useParPool,[],[],1);
assert( res(count)==saved.res(count) | (isnan(res(count)) & isnan(saved.res(count)) ))

%% Test 2, specified 2 [include bl, restrict]
count = 2;
i = 1;
thisClus = selClus(i);
theseSpikeTimes = spikeTimesAgg{thisClus};

useMaxDur = [-0.1 1];
restrictNeg = true; %estinate cannot be < 0
rng(1,'twister');
res(count) = latenzy(theseSpikeTimes,eventTimes,useMaxDur,...
    [],[],[],[],useParPool,[],restrictNeg);
assert( res(count)==saved.res(count) | (isnan(res(count)) & isnan(saved.res(count)) ))

%% Test 3, specified 4 [different jitter size, more resamples, use quantiles]
count = 3;
i = 1;
thisClus = selClus(i);
theseSpikeTimes = spikeTimesAgg{thisClus};

useMaxDur = 1;
resampNum = 500;
jitterSize = 4; %default is 2
useDirectQuant = true;
rng(1,'twister');
res(count) = latenzy(theseSpikeTimes,eventTimes,useMaxDur,resampNum,...
    jitterSize,[],[],useParPool,[],useDirectQuant);
assert( res(count)==saved.res(count) | (isnan(res(count)) & isnan(saved.res(count)) ))

%% Test 4, defaults
count = 4;
i = 2;
thisClus = selClus(i);
theseSpikeTimes = spikeTimesAgg{thisClus};

useMaxDur = 1;
rng(1,'twister');
res(count) = latenzy(theseSpikeTimes,eventTimes,useMaxDur,[],[],[],[],useParPool,[],[]);
assert( res(count)==saved.res(count) | (isnan(res(count)) & isnan(saved.res(count)) ))

%% Test 5, specified 2 [include bl, restrict]
count = 5;
i = 2;
thisClus = selClus(i);
theseSpikeTimes = spikeTimesAgg{thisClus};

useMaxDur = [-0.1 1];
restrictNeg = true; %estinate cannot be < 0
rng(1,'twister');
res(count) = latenzy(theseSpikeTimes,eventTimes,useMaxDur,...
    [],[],[],[],useParPool,[],restrictNeg);
assert( res(count)==saved.res(count) | (isnan(res(count)) & isnan(saved.res(count)) ))

%% Test 6, specified 4 [different jitter size, more resamples, use quantiles]
count = 6;
i = 2;
thisClus = selClus(i);
theseSpikeTimes = spikeTimesAgg{thisClus};

useMaxDur = 1;
resampNum = 500;
jitterSize = 4; %default is 2
useDirectQuant = true;
rng(1,'twister');
res(count) = latenzy(theseSpikeTimes,eventTimes,useMaxDur,resampNum,...
    jitterSize,[],[],useParPool,[],useDirectQuant);
assert( res(count)==saved.res(count) | (isnan(res(count)) & isnan(saved.res(count)) ))


%% Test 7, specified 4 [different jitter size, more resamples, use quantiles]
count = 7;
i = 3;
thisClus = selClus(i);
theseSpikeTimes = spikeTimesAgg{thisClus};

useMaxDur = 1;
resampNum = 500;
jitterSize = 4; %default is 2
useDirectQuant = true;
rng(1,'twister');
res(count) = latenzy(theseSpikeTimes,eventTimes,useMaxDur,resampNum,...
    jitterSize,[],[],useParPool,[],useDirectQuant);
assert( res(count)==saved.res(count) | (isnan(res(count)) & isnan(saved.res(count)) ))



