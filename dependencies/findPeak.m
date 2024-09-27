function [peakVal,peakTime] = findPeak(data,time,restrictRange,switchPosNeg)
% return highest peak (positive/negative) and corresponding timestamp, syntax:
%   [peakVal,peakTime] = findPeak(data,time,restrictRange,switchPosNeg)
%   input:
%   - data [N x 1]: values
%   - time [N x 1]: timestamps corresponding to data (default: 1:N)
%   - restrictRange [N x 2]: restrict peak to [restrictRange(1) restrictRange(2)] (default: [-inf inf])
%   - switchPosNeg: only detect (1) postive peak,(2) negative peak/trough, or (0) largest of the two (default: 0)
%
%   output:
%   - peakVal: peak value
%   - peakTime: time of peak
%
% history:
% 28 Aug 2023
%   - created by Robin Haak
% 30 Aug 2023
%   - changed from max(abs()) to peak detection
% 26 Sep 2023
%   - code clean-up
% 7 Feb 2024
%   - changed code to handle selection of the largest peak when one of the peaks vals is empty 
% 3 Mar 2024
%   - slight change in peak detection

%% prep
peakVal = nan;
peakTime = nan;

%check inputs
if numel(data) < 3
    return
end
if ~exist('time','var') || isempty(time)
    time = 1:numel(data);
end
if ~exist('restrictRange','var') || isempty(restrictRange)
    restrictRange = [-inf inf];
end
if ~exist('switchPosNeg','var') || isempty(switchPosNeg)
    switchPosNeg = 0;
end

%% find peaks
%postive
[valsPos,locsPos] = findpeaks(data);
indRem = time(locsPos)<restrictRange(1) | time(locsPos)>restrictRange(end);
valsPos(indRem) = [];
locsPos(indRem) = [];
[maxPos,idxMaxPos] = max(valsPos);
locMaxPos = locsPos(idxMaxPos);

%negative
[valsNeg,locsNeg] = findpeaks(-data);
indRem = time(locsNeg)<restrictRange(1) | time(locsNeg)>restrictRange(end);
valsNeg(indRem) = [];
locsNeg(indRem) = [];
[maxNeg,idxMaxNeg] = max(valsNeg);
locMaxNeg = locsNeg(idxMaxNeg);

%% select peak
loc = [];
if switchPosNeg == 1 && ~isempty(maxPos)
    %positive
    peakVal = maxPos;
    loc = locMaxPos;
elseif switchPosNeg == 2 && ~isempty(maxPos)
    %negative
    peakVal = -maxNeg;
    loc = locMaxNeg;
elseif switchPosNeg == 0
    %largest
    if ~isempty(maxPos) && ~isempty(maxNeg)
        if (maxPos > maxNeg)
            peakVal = maxPos;
            loc = locMaxPos;
        elseif maxPos < maxNeg
            peakVal = -maxNeg;
            loc = locMaxNeg;
        end
    elseif ~isempty(maxPos) && isempty(maxNeg)
        peakVal = maxPos;
        loc = locMaxPos;
    elseif isempty(maxPos) && ~isempty(maxNeg)
        peakVal = -maxNeg;
        loc = locMaxNeg;
    end
end

%return peak time
if ~isempty(loc)
    peakTime = time(loc);
end
