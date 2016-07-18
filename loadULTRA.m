function d = loadULTRAAnalysis(data_file_search_string,calibration_file_name)
d = construct2dPP;

files = dir(data_file_search_string);

temp = load(files.name);
calib = load(calibration_file_name);

d.w1 = fliplr(temp(1,1:end));
d.R  = fliplr(temp(2:end,:));
d.w3 = transpose(calib(1:end,2));
d.comment = sprintf('file %s calibration %s',files.name,calibration_file_name);

