function [sod, tod, fod, m , dmu_dq] = readGaussianAnharmonicCalculations(logname)
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
nModesString = '2nd derivatives larger than';
%quadraticString = ':      QUADRATIC FORCE CONSTANTS IN NORMAL MODES       :';
quadraticString = 'QUADRATIC FORCE CONSTANTS IN NORMAL MODES';
%cubicString = ':        CUBIC FORCE CONSTANTS IN NORMAL MODES         :';
cubicString = 'CUBIC FORCE CONSTANTS IN NORMAL MODES';
cubicEnumerateString = '3rd derivatives larger than';
%quarticString = ':       QUARTIC FORCE CONSTANTS IN NORMAL MODES        :';
quarticString = 'QUARTIC FORCE CONSTANTS IN NORMAL MODES';
quarticEnumerateString = '4th derivatives larger than';

indexOfSecondOrderDerivatives = 3;
indexOfThirdOrderDerivatives = 4;
indexOfFourthOrderDerivatives = 5;

quadraticHeaderLength = 8;
cubicHeaderLength = 8;
quarticHeaderLength = 8;

% search for the number of modes
hit = strfind(C{1}, nModesString);
nModesStart = find(~cellfun('isempty', hit));
thisline = C{1}(nModesStart);
thisline = thisline{1};
nmodes = sscanf(thisline,['%i ' nModesString '.*']);
if isempty(nmodes)
    tokens = regexp(thisline,'Num. of 2nd derivatives larger than .*:\s*(\d+) over.*','tokens');
    nmodes = str2num(tokens{1}{1});
end

% initialize second derivative vector
sod = zeros(1,nmodes);

% find the start of the quadratic section
hit = strfind(C{1}, quadraticString);
quadraticStart = find(~cellfun('isempty', hit)) + quadraticHeaderLength;

% looks like the quadratic constants are assumed to be diagonal (normal
% modes)
for ii  = 1:nmodes
    thisline = C{1}(quadraticStart+ii);
    thisline = thisline{1};
    stuff = sscanf(thisline,'%i %i %f %f %f');
    sod(ii) = stuff(indexOfSecondOrderDerivatives);
end


% lets do the third order now by reading the line, looking to see if it is
% empty or not, if it is empty break, if not process the line taking the
% coordinates from the first three numbers, assigning that value and looing

%find the start of the Cubic anharmonicity section
hit = strfind(C{1}, cubicString);
cubicStart = find(~cellfun('isempty', hit)) + cubicHeaderLength;

%find how many derivatives are listed in this section
hit = strfind(C{1},cubicEnumerateString);
thisline =  C{1}(find(~cellfun('isempty', hit)));
thisline = thisline{1};
n_tod = sscanf(thisline,['%i ' cubicEnumerateString '.*']);

%initialize tod
tod = zeros(n_tod,4);

modesDone = false;
count = 0;
while ~modesDone
    count = count+1;
    thisline = C{1}(cubicStart+count);
    thisline = thisline{1};
    if isempty(thisline)
        modesDone = true;
        break
    end
    stuff = sscanf(thisline,'%i %i %i %f %f %f');
%     ii = stuff(1);
%     jj = stuff(2);
%     kk = stuff(3);
%    tod(ii,jj,kk) = stuff(indexOfThirdOrderDerivatives);
    tod(count,1:4) = [stuff(1:3)' stuff(indexOfThirdOrderDerivatives)];
end

% lets do the fourth order now by reading the line, looking to see if it is
% empty or not, if it is empty break, if not process the line taking the
% coordinates from the first three numbers, assigning that value and looing

hit = strfind(C{1}, quarticString);
quarticStart = find(~cellfun('isempty', hit)) + quarticHeaderLength;

%find how many derivatives are listed in this section
hit = strfind(C{1},quarticEnumerateString);
thisline =  C{1}(find(~cellfun('isempty', hit)));
thisline = thisline{1};
n_fod = sscanf(thisline,['%i ' quarticEnumerateString '.*']);

%initialize fourth order derivatives matrix
fod = zeros(n_fod,5);

modesDone = false;
count = 0;
while ~modesDone
    count = count+1;
    thisline = C{1}(quarticStart+count);
    thisline = thisline{1};
    if isempty(thisline)
        modesDone = true;
        break
    end
    stuff = sscanf(thisline,'%i %i %i %i %f %f %f');
%     ii = stuff(1);
%     jj = stuff(2);
%     kk = stuff(3);
%     ll = stuff(4);
%     fod(ii,jj,kk,ll) = stuff(indexOfFourthOrderDerivatives);
    fod(count,1:5) = [stuff(1:4)' stuff(indexOfFourthOrderDerivatives)];
end

% search for the masses

mask = ~cellfun(@isempty, regexp(C{1}, 'Red. masses --'));
some_lines = C{1}(mask);
float_pat = '[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)';
matches = regexp(some_lines,float_pat,'match');
count = 0; 
areWeDone = 0;
m = zeros(1,nmodes);
for ii = 1:length(matches)
    if ~areWeDone
        for jj = 1:length(matches{ii})
            count = count+1;
            if count<=nmodes
                m(count) = str2double(matches{ii}{jj});
            else
                areWeDone = 1;
            end
        end
    end
end
m = fliplr(m);  % modes seem to be in the opposite order!

% search for the dipoles
%start_string = '## AFTER VARIATIONAL CORRECTION ##';
start_string = '## DEPERTURBED TRANSITIONS MOMENTS ##';
%stopString = 'Electric dipole : Overtones';
mask = ~cellfun(@isempty, regexp(C{1}, start_string));
ind = find(mask);
extra_lines = 6;
some_lines = C{1}(ind:ind+extra_lines+nmodes);
exponential_pat = '[-+]?[0-9]*\.?[0-9]+([eEdD][-+]?[0-9]+)';
mask = ~cellfun(@isempty,regexp(some_lines,exponential_pat,'match'));
some_lines = some_lines(mask);
matches = regexp(some_lines,exponential_pat,'match');
count = 0; 
areWeDone = 0;
dmu_dq = zeros(nmodes,3);
for ii = 1:length(matches)
    if ~areWeDone
        count = count+1;
        if count<=nmodes
            dmu_dq(count,1:3) = [str2double(replace(matches{ii}{1},'D','e')),...
                str2double(replace(matches{ii}{2},'D','e')),...
                str2double(replace(matches{ii}{3},'D','e'))];
        else
            areWeDone = 1;
        end
    end
end