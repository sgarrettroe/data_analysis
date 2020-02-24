function s = process2d(s, number, varargin)
% process2d will process the data and make a plot. 
% This should be enough for 90% of all needs 
%
% INPUT:
% - s (struct): the 2d-struct with all data
%
% - number (number): will be used as a number for the figure
%
% OPTIONAL INPUT, for details, see the functions they refer to.
% - apodization (name): will apply apodization in absorptive2d. Default is
%   'none'.
% - apod_nr ([a b]): factors for some of the apodization functions in
%   absorptive2d. Default is [10 10], which works for 'rbOnes' and
%   'rbGauss'.
% - debug (BOOL): will write some stuff, handy for debugging :)
% - fft_type (name): determines the type of fft in absorptive2d. Default is
%   'sgrsfft'.
% - n_contours (number): number of contours plotted by rb2dPlot. Default is
%   12.
% - zlimit (number): how much of the z-axis is plotted by rb2dPlot. Default
%   is 0 (all). 
% - no_units (BOOL): will print the indices instead of frequencies in
%   rb2dPlot. Handy for finding broken pixels etc. Default is false.


apodization = 'none';
flag_debug = false;
fft_type = 'sgrsfft';
n_contours = 12;
zlimit = 0;
no_units = false;
apod_nr = [-0.5 3];
line_width = 1;
xlim = [0 0];
ylim = [0 -1];

while length(varargin) >= 2
  arg = varargin{1};
  val = varargin{2};
  
  switch lower(arg)
    case 'apodization'
      apodization = val;
    case 'apod_nr'
      apod_nr = val;
    case 'xlim'
      xlim = val;
    case 'ylim'
      ylim = val;
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
    case 'line_width'
      line_width = val;
    otherwise
      error(['process2d: unknown option ', arg]);
  end 
  varargin = varargin(3:end);
end

if flag_debug; disp(['process2d: flag_debug:' int2str(flag_debug)]); end

s = absorptive2d(s,...
  'phase', s.phase,...
  'zeropad', s.zeropad,...
  'range', [s.freq(1), s.freq(end)],...
  'fft_type', fft_type,...
  'apodization', apodization,...
  'apod_numbers', apod_nr,...
  'debug', flag_debug,...
  'plot', flag_debug);

title_string = ({strcat(int2str(s.date),...
  ': ', s.basename,...
  ', T2: ', num2str(s.t2),...
  ', ts:', num2str(s.time_stamp),...
  ', scans x shots:', num2str(s.n_scan), ' x ', num2str(s.n_shots)),...
  s.comment}); 

figure(number),clf
rb2dPlot(s,... 
  'n_contours', n_contours,...
  'zlimit', zlimit,...
  'title', title_string,...
  'no_units', no_units,...
  'debug', flag_debug,...
  'line_width', line_width,...
  'xlim', xlim,...
  'ylim', ylim);


