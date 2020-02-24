function [L,R,T_degC,T_K] = loadULTRA(data_file_search_string,calibration_file_name)
n_pixels_per_array = 128;

path_name = fileparts(data_file_search_string);
files = dir(data_file_search_string);


calib = load(calibration_file_name);
w3 = calib(:,2)';

T_degC = zeros(length(files));
for ii = 1:length(files)
    L(ii) = construct2dPP;
    R(ii) = construct2dPP;
    
    str = fullfile(path_name,files(ii).name);
    temp = load(str);
    
    t = temp(1,2:end); %take first row and drop first dummy point
    temp = temp(2:end,2:end);%trim header row and column

    L(ii).w1 = [];
    L(ii).PP  = temp(1:n_pixels_per_array,:);
    L(ii).w3 = w3;
    L(ii).freq = w3;
    L(ii).comment = sprintf('file %s calibration %s',files(ii).name,calibration_file_name);
    L(ii).time = t;
    L(ii).time_units = 'ps';
    
    R(ii).w1 = [];
    R(ii).PP  = temp(n_pixels_per_array + (1:n_pixels_per_array),:);
    % last point is... different... throw it out for now
    R(ii).PP = R(ii).PP(1:end-1,:);
    R(ii).w3 = w3(1:end-1);
    R(ii).freq = w3(1:end-1);
    R(ii).comment = sprintf('file %s calibration %s',files(ii).name,calibration_file_name);
    R(ii).time = t;
    R(ii).time_units = 'ps';
    
    %try to infer t2 time from the file name
    t2 = [];
    uh = regexp(str,'([0-9]+) C\D*([0-9]+)fs.csv','tokens');
    if length(uh{1})==2
        T_degC(ii) = sscanf(uh{1}{1},'%i');
        t2     = sscanf(uh{1}{2},'%i');
    end
    L(ii).t2 = t2(1); %take only the first match?
    R(ii).t2 = t2(1);
    
    t0_bin = find(t==0);
    L(ii).t0_bin = t0_bin;
    R(ii).t0_bin = t0_bin;
    
    L(ii).zeropad = 4*length(t);
    L(ii) = absorptive2dPP(L(ii));
    R(ii) = absorptive2dPP(R(ii));
end
T_K = T_degC + 273;