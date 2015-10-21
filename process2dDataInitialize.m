if ispc
    username = getenv('Username');
elseif ismac
    username = getenv('User');
end

switch username
    case 'INFRARED' % 2D Lab computer
        datadirbasename = 'c:\data\';
        privatedir = 'C:\Users\INFRARED\Documents\GitHub\Mahou\spectrometer\private';
    case 'sgr' % Sean's computer
        datadirbasename = '/Volumes/data/';
        privatedir = '/Users/sgr/Lab/Mahou/spectrometer/private';
    case 'Tom Brinzer' % Tom on the Dell T3500
        datadirbasename = '\\Infrared-HP\data'; % path to the lab computer. May need to change if doing local work
        privatedir = 'C:\Users\Tom Brinzer\Documents\GitHub\Mahou\spectrometer\private';
    case 'Thomas Brinzer' % Tom on the updated laptop
        datadirbasename = '\\Infrared-HP\data';
        privatedir = 'C:\Users\Thomas Brinzer\Documents\GitHub';
        localdirbasename = strcat('c:\Users\',username,'\Documents\LocalData\');
    case 'Samrat' % Samrat on the Dell T3500
        datadirbasename = '\\Infrared-HP\data'; % path to the lab computer. May need to change if doing local work
        privatedir = 'C:\Users\Samrat\Documents\GitHub\Mahou\spectrometer\private';
    otherwise
        error('Unknown user. Unable to build data directory structure. Please add your info.');
end

% Build a full file directory from our input
datadir = fullfile(datadirbasename,datestring);
localdir = fullfile(localdirbasename,datestring);
tempdir = fullfile(datadirbasename,datestring,'temp');