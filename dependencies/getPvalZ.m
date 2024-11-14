function [pVals,Z] = getPvalZ(maxD,maxRandD,useDirectQuant)
%compute (z-scored) p-value, syntax:
%   [pVal,Z] = getPvalZ(maxD,maxRandD,useDirectQuant)
%
%based on code by Jorrit Montijn

%check inputs
if ~exist('useDirectQuant','var') || isempty(useDirectQuant)
    useDirectQuant = false;
end

%% calculate significance
%find highest peak and retrieve value
maxRandD = sort(unique(maxRandD));
randMu = mean(maxRandD);
randVar = var(maxRandD);

if useDirectQuant
    %calculate statistical significance using empirical quantiles
    pVals = nan(size(maxD));
    for i=1:numel(maxD)
        if maxD<min(maxRandD) || isnan(maxD) || isempty(maxRandD)
            value = 0;
        elseif maxD>max(maxRandD) || isinf(maxD) || numel(maxRandD)<3
            value = numel(maxRandD);
        else
            value = interp1(maxRandD,1:numel(maxRandD),maxD);
        end
        pVals(i) = 1-(value/(1+numel(maxRandD)));
    end

    %transform to output z-score
    Z = -norminv(pVals/2);
else
    %calculate statistical significance using Gumbel distribution
    [pVals,Z] = getGumbel(randMu,randVar,maxD);
end
end

