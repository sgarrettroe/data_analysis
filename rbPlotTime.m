function varargout = rbPlotTime(s, varargin)
%rbPlotTime plots the time domain spectra
%
%input is swith
%-data: R(1) or R(2)
%-pixel(opt): default is pixel 16. Pixel > 0 gives a trace of 1 pixel.
%pixel = 0 gives a 2d-plot
%
%output is a 2dplot 
%
%RB, 20110505: started function

pixel = 16;
xlim = [0 0];
ylim = [0 0];
while length(varargin)>=2
  arg = varargin{1};
  val = varargin{2};
  switch lower(arg)
    case 'plot'
      if strcmp(val, 'R1time')
        data = s.R1;
        axis_t = s.time;
        axis_w = s.freq;
      elseif strcmp(val, 'R2time')
        data = flipud(s.R2);
        axis_t = -flipud(s.time);
        axis_w = s.freq;
      elseif strcmp(val, 'noise')
        data = [flipud(s.R2_noise); s.R1_noise];
        axis_t = [-flipud(s.time); s.time];
        axis_w = s.freq; 
      elseif strcmp(val, 'time')
        % all time
        data = [flipud(s.R2); s.R1];
        axis_t = [-flipud(s.time); s.time];
        axis_w = s.freq;
      end
    case 'xlim'
      xlim = val;
    case 'ylim'
      ylim = val;
        
    case 'pixel'
      pixel = val;
      %see if pixel is in range
      shape = size(s.freq);
      if pixel < 0 && pixel > shape(2)
        error(['rbPlotTime: pixel is out of range ', arg])
      end
    otherwise
      error(['rbPlotTime: unknown option ',arg])
  end
  varargin = varargin(3:end);
end

%linear time trace
if pixel > 0
  plot(axis_t, data(:,1));
end

if pixel == 0
  rb2dPlot(axis_w, axis_t, data, 'xlim', xlim, 'ylim', ylim);
end



