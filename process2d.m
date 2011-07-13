function s = process2d(s, mess_date, number, varargin)
% process2d will process the data and make a plot. 
% This should be enough for 90% of all needs 


apodization = 'none';
flag_debug = false;
fft_type = 'sgrsfft';
n_contours = 12;
zlimit = 0;
no_units = false;
apod_nr = [10 10];

while length(varargin) >= 2
  arg = varargin{1};
  val = varargin{2};
  
  switch lower(arg)
    case 'apodization'
      apodization = val;
    case 'apod_nr'
      apod_nr = val;
    case 'debug'
      flag_debug = val;
    case 'fft_type'
      fft_type = val;
    case 'n_contours'
      n_contours = val;
    case 'zlimit'
      zlimit = val;
    case 'no_units'
      no_units = val;
    otherwise
      error(['process2d: unknown option ', arg]);
  end 
  varargin = varargin(3:end);
end

s = absorptive2d(s,...
  'phase', s.phase,...
  'zeropad', s.zeropad,...
  'range', [s.freq(1), s.freq(end)],...
  'fft_type', fft_type,...
  'apodization', apodization,...
  'apod_numbers', apod_nr,...
  'debug', flag_debug,...
  'plot', false);

title_string = ({strcat(num2str(mess_date),...
  ': ', s.basename,...
  ', T2: ', num2str(s.t2),...
  ', ts:', num2str(s.time_stamp),...
  ', scans x shots:', num2str(s.n_scan), ' x ', num2str(s.n_shots)),...
  s.comment}); 

figure(number),clf
rb2dPlot(s, 'n_contours', n_contours, 'zlimit', zlimit, 'title', title_string, 'no_units', no_units);


