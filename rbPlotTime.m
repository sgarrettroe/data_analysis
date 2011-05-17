function varargout = rbPlotTime(data, varargin)
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
while length(varargin)>=2
  arg = varargin{1};
  val = varargin{2};
  switch lower(arg)
    case 'pixel'
      pixel = val;
      %see if pixel is in range
      shape = size(data);
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
  plot(data(:,pixel))
end

%2dplot
if pixel == 0
  shape = size(data);
  x = linspace(0, shape(1), shape(1));
  y = linspace(1, shape(2), shape(2));
  
  my2dPlot(x,y, transpose(data),'n_contours', 12)
end