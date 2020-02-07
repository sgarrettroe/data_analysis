function [sod, tod, fod,m,dmu_dq] = readQChemAnharmonicCalculations(logname)
% read second, third, and fourth order derivatives from a Gaussian
% anharmonic vibrational analysis log file.
sod = [];
tod = [];
fod = [];
m = [];
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


% initialize second derivative vector
sod = zeros(1,nmodes);
m = zeros(1,nmodes);
dmu_dq = zeros(nmodes,3);

% find the start of the quadratic section
hits = strfind(C{1}, quadraticString);
wlines = find(~cellfun('isempty', hits));
hits = strfind(C{1}, massString);
mlines = find(~cellfun('isempty', hits));
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
    %extract transition dipoles
    thisline = C{1}(tlines(ii));
    thisline = thisline{1};
    tstuff = sscanf(thisline,[dmu_dqString '%f %f %f %f %f %f %f %f %f']);
    
    for jj = 1:length(wstuff)
        count = count+1;
        sod(count) = wstuff(jj);
        m(count) = mstuff(jj);
        dmu_dq(count,1:3) = tstuff([3*(jj-1)+1, 3*(jj-1)+2, 3*(jj-1)+3]);
    end
end

% search for third-order anharmonicities, strings like, e.g.
% Eta[1,1,1]=-63.9. Capture the indices and the value as tokens.
tokens = regexp(C{1},'Eta\[(\d+),(\d+),(\d+)\]=(-*\w*.\w*),','tokens');
hits = find(~cellfun('isempty',tokens));
tokens = tokens(hits);
n_tod = length(tokens);

%initialize tod
tod = zeros(n_tod,4);
for ii = 1:n_tod
    tod(ii,1) = sscanf(tokens{ii}{1}{1},'%i');
    tod(ii,2) = sscanf(tokens{ii}{1}{2},'%i');
    tod(ii,3) = sscanf(tokens{ii}{1}{3},'%i');
    tod(ii,4) = sscanf(tokens{ii}{1}{4},'%f');
end

% lets do the fourth order now 
tokens = regexp(C{1},'Eta\[(\d+),(\d+),(\d+),(\d+)\]=(-*\w*.\w*),','tokens');
hits = find(~cellfun('isempty',tokens));
tokens = tokens(hits);
n_fod = length(tokens);

%initialize tod
fod = zeros(n_fod,5);
for ii = 1:n_fod
    fod(ii,1) = sscanf(tokens{ii}{1}{1},'%i');
    fod(ii,2) = sscanf(tokens{ii}{1}{2},'%i');
    fod(ii,3) = sscanf(tokens{ii}{1}{3},'%i');
    fod(ii,4) = sscanf(tokens{ii}{1}{4},'%i');
    fod(ii,5) = sscanf(tokens{ii}{1}{5},'%f');
end
