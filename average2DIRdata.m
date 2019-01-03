function dataOut = average2DIRdata(dataIn,varargin)
% This code averages multiple 2D-IR spectra at a single time (t2) point. It
% does this by cropping the non-negative time elements in each spectrum (at
% a given t2) to the same length, then averaging them, and outputting a new
% data object.
%
% It is designed to work on a MatLab data object of the type we use for our
% 2D-IR experiments.
% 
% It will also not return a correct result if the same t1 spacing is not
% used. I suspect that it will not work at all, but you might get unlucky
% and have it work.
%
% You can currently use this code on 2D-IR spectra with different values of
% w3 (in terms of the angle). Thus, it is very dangerous to use on data
% from different days, unless you have ensured that the w3 values are the
% same. Similarly, you will get very wrong results if you use this with
% data where portions of it have been calibrated in w3 and portions have
% not. That being said, none of this has ever come up in day-to-day
% operation for me -- I just want to make sure that it's noted here that we
% don't check for w3 similarity currently.
%
% TODO: Check for w3 values.

n_zp = 1;
phase_deg = 0;
apod = 'none';
filter = false;

while length(varargin)>=2 %using a named pair
    arg = varargin{1};
    val = varargin{2};
    switch lower(arg)
        case 'n_zp'
            n_zp = val;
        case 'delta_phase'
            phase_deg = val;
        case 'apodization'
            apod = val;
        case 'filter'
            filter = val;
        otherwise
            warning(['average2DIRdata: unknown option ',arg])
    end
    varargin = varargin(3:end);
end

nonNegTime = zeros(size(dataIn)); % initialize the vector of non-negative times

% for ease of use later, create vectors containing the t2 times (in fs) and
% the t0_bins of each element in dataIn
t2s = [dataIn.t2];
t0_bin = [dataIn.t0_bin];

% Populate a vector containing the length of the PP data (in t1) from t0 to
% the end of the vector (the non-negative times)
for ii = 1:length(dataIn)
    [~,nonNegTime(ii)] = size(dataIn(ii).PP(:,t0_bin(ii):end));
end

shortestTimeLength = min(nonNegTime); % find the shortest non-negative time length

% crop all of the pump-probe data to be the same (non-zero) length
for ii = 1:length(dataIn)
    dataIn(ii).time = dataIn(ii).time(1:t0_bin(ii) + shortestTimeLength - 1);
    dataIn(ii).PP = dataIn(ii).PP(:,1:t0_bin(ii) + shortestTimeLength - 1);
end

% determine the unique population time points, and the index vector that
% maps their new positions to their old (non-unique) positions
[uniquet2s,~,ic] = unique(t2s);

for ii = 1:length(uniquet2s)
    index = ic==ii; % take only the elements that match our current t2
    temp = dataIn(index); % make a smaller version of data for processing
    zeropad = temp(1).zeropad;
    nscans = zeros(size(temp));
    for jj = 1:length(temp)
        temp(jj) = absorptive2dPP(temp(jj),'zeropad',n_zp*zeropad,...
            'phase',temp(jj).phase + phase_deg,'filter',filter,'apodization',apod); % re-run absorptive 2d PP to get correct w1 spacing
        temp(jj).zeropad = n_zp*zeropad;
        nscans(jj) = temp(jj).PARAMS.nScans;
    end
    [m,n] = size(temp(1).R); % figure out the new size of R
% figure out your average R -- we turn elements of the temp structure into
% a 3-D matrix, and then average along the third dimension
    R = mean(reshape([temp.R],m,n,length(temp)),3); 
    
    dataOut(ii) = temp(1);
    dataOut(ii).R = R;
    dataOut(ii).phase = [];
    dataOut(ii).PP = [];
    dataOut(ii).t2 = uniquet2s(ii);
    dataOut(ii).t0_bin = [];
    dataOut(ii).igram = [];
    dataOut(ii).PARAMS.nScans = sum(nscans);
    disp(ii);
end