function s = subtraction2d(plus, minus)
% subtraction2d subtracts two 2d structs and makes a new one
% it will check if certain parameters are the same (axes, ...) and give an
% error if they aren't
% it will check other parameters if they are the same (t2, ...) and give a
% warning if they aren't
% it will copy some parameters if they are the same (date, phase, ...)
% it will just copy some other parameters from "plus"


s = construct2d;
if nargin == 0
  return
end

% things that are different 
s.R = plus.R - minus.R;
s.basename = [plus.basename ' minus ' minus.basename];
s.comment = 'Subtraction of two spectra';

s.time_stamp = 0;

% things that don't have to be copied, if they are not the same
if plus.phase == minus.phase
  s.phase = plus.phase;
end

if plus.date == minus.date
  s.date = plus.date;
end

if plus.n_shots == minus.n_shots
  s.n_shots = plus.n_shots;
end

if plus.n_scan == minus.n_scan
  s.n_scan = plus.n_scan;
end

% things that have to be the same -> gives an error
if plus.w1 ~= minus.w1
  error(['The w1 axes of the two are not the same.'])
else 
  s.w1 = plus.w1;
end

if plus.w3 ~= minus.w3
  error(['The w3 of the two are not the same.'])
else 
  s.w3 = plus.w3;
end

if plus.centerfreq ~= minus.centerfreq
  error(['The center frequency of the two are not the same.'])
else 
  s.centerfreq = plus.centerfreq;
end

if plus.undersampling ~= minus.undersampling
  error(['The center frequency of the two are not the same.'])
else 
  s.undersampling = plus.undersampling;
end

if plus.resolution ~= minus.resolution
  error(['The resolution of the two are not the same.'])
else 
  s.resolution = plus.resolution;
end

if plus.zeropad ~= minus.zeropad
  error(['The zeropadding of the two are not the same.'])
else 
  s.zeropad = plus.zeropad;
end

% things that should be the same -> gives a warning
if plus.freq ~= minus.freq
  warning(['Frequencies of the two are not the same.'])
else
  s.freq = plus.freq;
end

if plus.time ~= minus.time
  warning(['The time steps of the two are not the same.'])
else 
  s.time = plus.time;
end

if plus.t2 ~= minus.t2
  warning(['The population time of the two are not the same.'])
else 
  s.t2 = plus.t2;
end

if plus.t3 ~= minus.t3
  warning(['The t3 of the two are not the same.'])
else 
  s.t3 = plus.t3;
end

if plus.spec_calib ~= minus.spec_calib
  warning(['The spectral calibration of the two are not the same.'])
else 
  s.spec_calib = plus.spec_calib;
end

% things that are unlikely to be different and I'm lazy
s.time_units = plus.time_units;
s.freq_units = plus.freq_units;
s.pump_probe = plus.pump_probe;
s.pump_probe_freq = plus.pump_probe_freq;





















