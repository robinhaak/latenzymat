[spikeTimes,eventTimes,onsetTimes,minLatency,useMaxDur,sponRates,peakRates,peakDurs,sustRates,jitterSize] = ...
    generateTrialSpikesPeak(50,50,0.1);

set(0,'DefaultFigureWindowStyle','docked')

for thisUnit = 1:numel(spikeTimes)
    theseSpikeTimes = spikeTimes{thisUnit};
     [respLatency,sLatenzy] = ...
         latenzy(theseSpikeTimes,eventTimes,[0 2],[],[],[],[],[],[],[],1);
     drawnow;
end



