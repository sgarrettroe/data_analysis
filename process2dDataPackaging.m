%% Initialize
process2dDataInitialize

%% Load data and build MatLab structure

 % 1xn vector so it plays nicely with the FOR loop later
ndirs = length(dirindex);
scanoffsetlist = cell(1,ndirs); % To keep track of our scan offset.
% Builds a 1 x ndirs cell array of empty matrices

clear data c2
scanoffset = 0; % So we can sew our data together later ... see ~ line 85

% Build of data array
for ii = 1:ndirs
    %choose directory
    thisdatadir = datadir{dirindex(ii)};
    thisscanlist_holder = thisscanlist{ii}; %curly braces pull the elements
    scanoffsetlist{ii} = scanoffset; %
    
    cd(thisdatadir) %move into directory
    
    % Can we pre-allocated space for the array? Probably not.
    for jj = thisscanlist_holder
        %Loads data from our *.mat file
        m = load(sprintf('%03i.mat',jj)); % Why the syntax?
        data(scanoffset + jj) = m.data;
        data(scanoffset + jj).w3 = spectrum_calibration(1).*data(scanoffset ...
            + jj).w3 + spectrum_calibration(2);
        data(scanoffset + jj).spec_calib = spectrum_calibration;
    end
    for jj = thisscanlist_holder
        data(scanoffset + jj).scan_number = scanoffset + jj;
        data(scanoffset + jj).datestring = datestring;
        data(scanoffset + jj).name = name;
        data(scanoffset + jj).temperature = temperature;
    end
    scanoffset = length(data); %Now I see it. This allows us to splice together 
    % different days into the same array, without loss of indexing. Think
    % of it as sewing ...
end
%% Remove empty elements from the data structured array, and sort by t2
empty_elems = arrayfun(@(s) all(structfun(@isempty,s)),data);
data2 = data(~empty_elems);
clear data
uh = [data2.t2];
[~,ind] = sort(uh);
data = data2(ind);
clear data2 ind uh empty_elems
%% Initial data visualization
% This allows us to check and see if the quality of our spectra is good,
% and if our calibration is correct.
for kk = 1:length(data);%thisscanlist_holder
    % display figure
    ind1 = find(data(kk).w1>range1(1) & data(kk).w1<range1(2));
    ind3 = find(data(kk).w3>range3(1) & data(kk).w3<range3(2));
    figure(2),clf,
    my2dPlot(data(kk).w1(ind1),data(kk).w3(ind3),data(kk).R(ind3,ind1),...
        'pumpprobe',false,'n_contours',20)
    % add 'zlimit' flag to my2dPlot if you need to zoom in
    pause
end

%% Save data file
beep

flag_save = input('Save output? (1 = "YES"; 0 = "NO"): ');
if flag_save
    cd(localdir{1})
    save(sprintf('%s-%02i_summary.mat',datestring{dirindex(1)},experiment{1}),'data', ...
        'datestring','datadirbasename','experiment','range1','range3','ind1','ind3')
end
flag_save = 0;