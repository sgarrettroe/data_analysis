function dataOut = average2DIRdata(dataIn)
% This code averages multiple 2D-IR spectra at a single time (t2) point. It
% does this by cropping the non-negative time elements in each spectrum (at
% a given t2) to the same length, then averaging them, and outputting a new
% data object.
%
% It is designed to work on a MatLab data object of the type we use for our
% 2D-IR experiments.
% 
% It will break horrendously if the same t1 spacing is not used ... so
% caveat emptor.


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
    for jj = 1:length(temp)
        temp(jj) = absorptive2dPP(temp(jj)); % re-run absorptive 2d PP to get correct w1 spacing
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
end