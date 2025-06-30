% testLatenzy2
%
% To run use runtests('testLatenzy2')
%
% 2025, Robin Haak, Alexander Heimel

clear; rng(1,'twister');

pth = fileparts(which('testLatenzy'));
load(fullfile(pth,'Topo2_20220126_AP.mat'));

idxIncl = ([sAP.sCluster.KilosortGood] | [sAP.sCluster.Contamination] < 0.1) & ...
    contains({sAP.sCluster.Area},'primary visual','IgnoreCase',true);

spikeTimesAgg = {sAP.sCluster(idxIncl).SpikeTimes};
numClus = numel(spikeTimesAgg);

sBlock = sAP.cellBlock{4}; %DG, 1s
eventTimes = sBlock.vecStimOnTime;

%create two conditionsgen
n = numel(eventTimes);
idx = my_randperm(n,n);
eventTimes1 = sort(eventTimes(idx(1:floor(n/2))));
eventTimes2 = sort(eventTimes(idx(floor(n/2)+1:end)));
%jitter = zeros(size(eventTimes2));
%numJitter = round(3/4 * numel(eventTimes2));
%jitterIdx = randperm(numel(eventTimes2), numJitter);
%jitter(jitterIdx) = (rand(numJitter,1)*2 - 1) * 2;  % ±2 sec jitter
%eventTimes2 = eventTimes2 + jitter;

useParPool = false;
selClus = [16 19 20];
numClus = length(selClus);

% move spikes from eventTimes1 to eventTimes2
diff_latency = 0.1;
for i=1:numClus
    thisClus = selClus(i);
    theseSpikeTimes = spikeTimesAgg{thisClus};

    for j=1:length(eventTimes1)
        ind = find(theseSpikeTimes>(eventTimes1(j)+diff_latency) & theseSpikeTimes<(eventTimes1(j)+diff_latency+0.1));
        if ~isempty(ind)
            % select a quarter of these spikes 
            n = length(ind);
            k = ceil(0.25*n);
            ind = ind(my_randperm(n,k));
            % and move them to the second set of events
            if ~isempty(ind)
                newSpikeTimes = theseSpikeTimes(ind)- eventTimes1(j) + eventTimes2(j);
                theseSpikeTimes(ind) = [];
                theseSpikeTimes = sort([theseSpikeTimes; newSpikeTimes]);
            end
        end
    end
    spikeTimesAgg{thisClus} = theseSpikeTimes;
    find(theseSpikeTimes>eventTimes1(1)+0.1 &theseSpikeTimes<eventTimes1(1)+0.2)
end
res = NaN(1,6); 
saved = load(fullfile(pth,'testResults2.mat'),'res');

%% Test 1 - default
count = 1;
i = 1;
thisClus = selClus(i);
theseSpikeTimes = spikeTimesAgg{thisClus};
useMaxDur = 1;
makePlots = 1;
rng(1,'twister');
res(count) = latenzy2(theseSpikeTimes,eventTimes1,theseSpikeTimes,eventTimes2,useMaxDur,...
    [],[],useParPool,[],[],makePlots);
%assert( res(count)==saved.res(count) | (isnan(res(count)) & isnan(saved.res(count)) ))

%% Test 2 - specified [include bl, restrict]
count = 2;
i = 1;
thisClus = selClus(i);
theseSpikeTimes = spikeTimesAgg{thisClus};
useMaxDur = [-0.1 1];
restrictNeg = true;
rng(1,'twister');
res(count) = latenzy2(theseSpikeTimes,eventTimes1,theseSpikeTimes,eventTimes2,useMaxDur,...
    [],[],useParPool,[],restrictNeg);
%assert( res(count)==saved.res(count) | (isnan(res(count)) & isnan(saved.res(count)) ))

%% Test 3 - specified [more resamps, use quantiles]
count = 3;
i = 1;
thisClus = selClus(i);
theseSpikeTimes = spikeTimesAgg{thisClus};
useMaxDur = 1;
resampNum = 500;
useDirectQuant = true;
rng(1,'twister');
res(count) = latenzy2(theseSpikeTimes,eventTimes1,theseSpikeTimes,eventTimes2,useMaxDur,...
    resampNum,[],useParPool,useDirectQuant,[]);
%assert( res(count)==saved.res(count) | (isnan(res(count)) & isnan(saved.res(count)) ))

%% Test 4 - unequal trial N
count = 4;
i = 1;
thisClus = selClus(i);
theseSpikeTimes = spikeTimesAgg{thisClus};
useMaxDur = 1;
trialIdx1 = 1:3:numel(eventTimes1);
trialIdx2 = 1:4:numel(eventTimes1);
rng(1,'twister');
res(count) = latenzy2(theseSpikeTimes,eventTimes1(trialIdx1),theseSpikeTimes,...
    eventTimes2(trialIdx2),useMaxDur,[],[],useParPool,[],[]);
%assert( res(count)==saved.res(count) | (isnan(res(count)) & isnan(saved.res(count)) ))

%% Test 5 - default
count = 5;
i = 2;
thisClus = selClus(i);
theseSpikeTimes = spikeTimesAgg{thisClus};
useMaxDur = 1;
makePlots = 0;
rng(1,'twister');
res(count) = latenzy2(theseSpikeTimes,eventTimes1,theseSpikeTimes,eventTimes2,useMaxDur,...
    [],[],useParPool,[],[],makePlots);
%assert( res(count)==saved.res(count) | (isnan(res(count)) & isnan(saved.res(count)) ))

%% Test 6 - default
count = 6;
i = 3;
thisClus = selClus(i);
theseSpikeTimes = spikeTimesAgg{thisClus};
useMaxDur = 1;
makePlots = 0;
rng(1,'twister');
res(count) = latenzy2(theseSpikeTimes,eventTimes1,theseSpikeTimes,eventTimes2,useMaxDur,...
    [],[],useParPool,[],[],makePlots);
%assert( res(count)==saved.res(count) | (isnan(res(count)) & isnan(saved.res(count)) ))


