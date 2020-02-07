function data_out = calibrate2DIRdata(data_in,calibration)
% calibrate2DIRdata - Calibrates the w3 axis of a 2D-IR spectrum.
% 'calibration' should be a row vector contain the slope and offset of the
% calibration. Currently, this is only set up for a linear calibration.
%   
%      data = calibration2DIRdata(data,[1 0]);
%

for ii = 1:length(data_in)
        data_in(ii).w3 = calibration(1).*data_in(ii).freq + calibration(2);
        data_in(ii).spec_calib = calibration;
end
data_out = data_in;