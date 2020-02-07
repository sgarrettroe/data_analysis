function [w,intensity,m,dmu_dq] = readQChemHarmonicCalculations(logname)
% read second, third, and fourth order derivatives from a Gaussian
% anharmonic vibrational analysis log file.
w = [];
m = [];
intensity = [];
dmu_dq = [];

if exist(logname,'file')
    fid = fopen(logname,'r');
else
    error(['couldn''t find ',logname,' in search path. Check working directory']);
end

C = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);

%// Search a specific string and find all rows containing matches
%quadraticString = ':      QUADRATIC FORCE CONSTANTS IN NORMAL MODES       :';
quadraticString = 'Frequency:';
massString = 'Red. Mass:';
intensityString = 'IR Intens:';
dmu_dqString = 'TransDip';

% search for the number of modes
nModesString = 'Mode:';
hit = strfind(C{1}, nModesString);
nModesStart = find(~cellfun('isempty', hit));
%uh = regexp(C{1}(nModesStart),'Mode:\s*(\d+\s*)+','tokens')
nmodes = 0;
for ii = 1:length(nModesStart) 
    thisline = C{1}(nModesStart(ii));
    thisline = thisline{1};
    uh = sscanf(thisline,'Mode:%i %i %i');
    if max(uh)>nmodes
        nmodes = max(uh);
    end
end


% initialize vectors
w = zeros(1,nmodes);
m = zeros(1,nmodes);
intensity = zeros(1,nmodes);
dmu_dq = zeros(nmodes,3);

% find the start of the quadratic section
hits = strfind(C{1}, quadraticString);
wlines = find(~cellfun('isempty', hits));
hits = strfind(C{1}, massString);
mlines = find(~cellfun('isempty', hits));
hits = strfind(C{1}, intensityString);
ilines = find(~cellfun('isempty', hits));
hits = strfind(C{1}, dmu_dqString);
tlines = find(~cellfun('isempty', hits));

% looks like the quadratic constants are assumed to be diagonal (normal
% modes)
count = 0;
for ii  = 1:ceil(nmodes/3)
    %extract frequencies
    thisline = C{1}(wlines(ii));
    thisline = thisline{1};
    wstuff = sscanf(thisline,[quadraticString '%f %f %f']);
    %extract reduced masses
    thisline = C{1}(mlines(ii));
    thisline = thisline{1};
    mstuff = sscanf(thisline,[massString '%f %f %f']);
    %extract intensities
    thisline = C{1}(ilines(ii));
    thisline = thisline{1};
    istuff = sscanf(thisline,[intensityString '%f %f %f']);
    %extract transition dipoles
    thisline = C{1}(tlines(ii));
    thisline = thisline{1};
    tstuff = sscanf(thisline,[dmu_dqString '%f %f %f %f %f %f %f %f %f']);
    
    for jj = 1:length(wstuff)
        count = count+1;
        w(count) = wstuff(jj);
        m(count) = mstuff(jj);
        intensity(count) = istuff(jj);
        dmu_dq(count,1:3) = tstuff([3*(jj-1)+1, 3*(jj-1)+2, 3*(jj-1)+3]);
    end
end
