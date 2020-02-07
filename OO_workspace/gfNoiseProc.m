function dataStruct2 = gfNoiseProc(dataStruct,basedatadir)
% For now, this is not a generalized search. It only works when you have
% your data file system set up to have folders by date (YYYY-MM-DD), with a
% /temp directory under them. It also only works on Windows machines (I
% think) because of the direction of the slashes. 
homedir = pwd;

for ii = 1:length(dataStruct)
    count = 0;
    % Gives us the set of temp files to go to
    run_number = dataStruct(ii).scan_number;

    % requires your data to have a datestring, but puts us into the correct
    % temp directory if your data is split across multiple days
    datestring = datestr(dataStruct(ii).datestring,29);
    datadir = fullfile(basedatadir,datestring,'\temp');
    cd(datadir)
    
    % defines the step for our scattering suppression. For now, this is
    % two, but could change with a different scattering suppression method.
    phase_steps = 2;
    
    % Total number of scans we need to process through.
    nScans = dataStruct(ii).PARAMS.nScans;
    calibration = dataStruct(ii).spec_calib;
    if nScans < 4
        error('Your spectrum has too few scans for this script to work.');
    end
    
    % Define a stating point for the response (R). We'll redefine this
    % repeatedly later.
    oldR = dataStruct(ii).R;
    
    % Initialize matrix where we'll accumulate R for each set of two scans.
    dR = zeros([size(oldR),floor(nScans/phase_steps)]); 
    
    % Get the range to crop our data from temp files over.
    w1_range = [dataStruct(ii).w1(1) dataStruct(ii).w1(end)];
    w3_range = [dataStruct(ii).w3(1) dataStruct(ii).w3(end)];
    
    % ESTIMATE THE NOISE
    %
    % Load each spectrum from the temp file, and from every intensity in
    % each temp file, calculate the standard deviation 
    %
    % The temp files are cumulative, so to get an individual spectrum you have to subtract
    % the previous one.
    % 
    % The scattering suppression (done by moving the population stage) is
    % only good for even numbered scans, and scattering can add to the
    % standard deviation. Therefore, we only take every other temp file, so
    % that we have only spectra with scattering suppressed.
    
    for iii = nScans - phase_steps:-phase_steps:1
        fprintf('%s %i; %s %i; %i %s %i\n','Run',run_number,' Scan',iii,count,'of',floor(nScans/phase_steps))
        count = count + 1;

        % we took care of CDing to the right directory earlier
        m = load(sprintf('%03i-%04i.mat',run_number,iii));
        m.data = calibrate2DIRdata(m.data,calibration);
        % Excise the portion of the spectrum that matches our input
        % structure
        cropped_temp = cropData(m.data,w1_range,w3_range);
        
        % Build a matrix of the back-calculated spectral intensities for
        % each set of two scans. We weight each spectrum by the scan number
        % to get the intensities, since it is an average (rather than sum)
        % of the preceding data.
        dR(:,:,count) = oldR.*(iii + phase_steps) - cropped_temp.R*iii;
        
        % Redefine our 'oldR' for the next cycle of the loop
        oldR = cropped_temp.R;
    end
    
    % calculate the standard deviation for each data point across all of
    % the scans we processed.

    dataStruct(ii).sd = std(dR,[],3);
    dataStruct(ii).rsd = dataStruct(ii).sd./(mean(abs(dR),3));
%     dataStruct(ii).vmr = dataStruct(ii).sd.^2../(mean(abs(dR),3)); 
%     
    dataStruct2 = dataStruct;
    disp('Done')
end
cd(homedir)