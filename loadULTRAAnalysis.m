function d = loadULTRAAnalysis(data_file_search_string,calibration_file_name)
d = construct2dPP;

path_name = fileparts(data_file_search_string);
files = dir(data_file_search_string);

for ii = 1:length(files)
    str = fullfile(path_name,files(ii).name);
    temp = load(str);
    calib = load(calibration_file_name);
    
    d(ii).w1 = fliplr(temp(1,1:end));
    d(ii).R  = fliplr(temp(2:end,:));
    d(ii).w3 = transpose(calib(1:end,2));
    d(ii).comment = sprintf('file %s calibration %s',files(ii).name,calibration_file_name);
    
    %try to infer t2 time from the file name
    t2 = sscanf(files(ii).name,'%i');
    d(ii).t2 = t2(1); %take only the first match?
end
