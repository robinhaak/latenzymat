%testLatenzy_20250528
%Robin Haak
clear; rng(1,'twister');

%% load data
try
    load(fullfile(pwd,'Topo2_20220126_AP.mat'));
catch ME
    disp(ME.message);
end

idxIncl = ([sAP.sCluster.KilosortGood] | [sAP.sCluster.Contamination] < 0.1) & ...
    contains({sAP.sCluster.Area},'primary visual','IgnoreCase',true);

spikeTimesAgg = {sAP.sCluster(idxIncl).SpikeTimes};
numClus = numel(spikeTimesAgg);

sBlock = sAP.cellBlock{4}; %DG, 1s
eventTimes = sBlock.vecStimOnTime;

%% run tests
useParPool = false;
res = nan(5,numClus); %n_tests x n_clus

for thisClus = 1:numClus

    fprintf('Running tests for cluster %d/%d\n',thisClus,numClus);
    theseSpikeTimes = spikeTimesAgg{thisClus};

    %default
    useMaxDur = 1;
    rng(1,'twister');
    res(1,thisClus) = latenzy(theseSpikeTimes,eventTimes,useMaxDur,[],[],[],[],useParPool,[],[]);

    %specified 1 [include bl]
    useMaxDur = [-0.1 1];
    restrictNeg = false;
    rng(1,'twister');
    res(2,thisClus) = latenzy(theseSpikeTimes,eventTimes,useMaxDur,...
        [],[],[],[],useParPool,[],restrictNeg);

    %specified 2 [include bl, restrict]
    useMaxDur = [-0.1 1];
    restrictNeg = true; %estinate cannot be < 0
    rng(1,'twister');
    res(3,thisClus) = latenzy(theseSpikeTimes,eventTimes,useMaxDur,...
        [],[],[],[],useParPool,[],restrictNeg);

    %specified 3 [more resamps, use quantiles]
    useMaxDur = 1;
    resampNum = 500;
    useDirectQuant = true;
    rng(1,'twister');
    res(4,thisClus) = latenzy(theseSpikeTimes,eventTimes,useMaxDur,resampNum,...
        [],[],useParPool,[],useDirectQuant,[]);

    %specified 4 [different jitter size]
    useMaxDur = 1;
    jitterSize = 4; %default is 2
    rng(1,'twister');
    res(4,thisClus) = latenzy(theseSpikeTimes,eventTimes,useMaxDur,resampNum,...
        jitterSize,[],[],useParPool,[],[]);
end



