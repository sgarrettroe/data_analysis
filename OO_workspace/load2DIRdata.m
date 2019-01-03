function output = load2DIRdata(directory,scanlist)
% struct = load2DIRdata(directory,scanlist)
%
% This function loads 2D-IR data from one or more directories into MatLab's
% stack. 'Directory' should be a string. 'Scanlist' should be a vector
% listing the scans you want to import. If your data is split across
% multiple directories, you can use cell arrays for 'directory' and
% 'scanlist'. An example of the syntax is:
%
%     data = load2DIRdata('\\Infrared-HP\data\1989-07-25',10:20);
%         OR
%     data = load2DIRdata({'\\Infrared-HP\data\1989-07-25', ...
%                        '\\Infrared-HP-data\1989-07-26'},{10:20,[1:2 5:8]}
%
% The data pulled using this command will be UNSORTED, NON-PHASE CORRECTED,
% and UNCALIBRATED. You'll probably need to fix this before doing any
% analysis on it.
old_dir = pwd;

if ~isa(directory,'cell')
    directory = {directory};
end
if ~isa(scanlist,'cell')
    scanlist = {scanlist};
end

ndirs = length(directory);

% Builds a 1 x ndirs cell array of empty matrices
scanoffsetlist = cell(1,ndirs); % To keep track of our scan offset.

% Initialize the SCANOFFSET variable that we use to  patch sets of data
% from multiple directories together
scanoffset = 0; 
totalscans = 0;

for ii = 1:ndirs
    totalscans = totalscans + length(scanlist{ii});
    datestring = cell(1,totalscans);
end


% Build data array
for ii = 1:ndirs
    %choose directory
    datadir = directory{ii};
    scanlist_holder = scanlist{ii}; %curly braces pull the elements
    scanoffsetlist{ii} = scanoffset; %
    
    cd(datadir) %move into directory
    
    % Can we pre-allocated space for the array? Probably not.
    for jj = scanlist_holder
        %Loads data from our *.mat file
        filename = sprintf('%03i.mat',jj);
        D = dir(filename);
        datestring{scanoffset + jj} = D.date;
        m = load(filename);
        output(scanoffset + jj) = m.data;
    end
    for jj = scanlist_holder
        output_(scanoffset + jj).scan_number = jj;
        output_(scanoffset + jj).datestring = datestring{scanoffset + jj};
%         data(scanoffset + jj).name = name;
%         data(scanoffset + jj).temperature = temperature;
    end
    scanoffset = length(output); %Now I see it. This allows us to splice together 
    % different days into the same array, without loss of indexing.
    cd(old_dir);
end
for ii = 1:length(output)
    output(ii).scan_number = output_(ii).scan_number;
    output(ii).datestring = output_(ii).datestring;
end