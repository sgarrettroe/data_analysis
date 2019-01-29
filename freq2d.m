function s = freq2d(s,varargin)
%freq2d set up the frequency axis with undersampling etc
%function s = freq2d(s,[zeropad])


flagSpectrometer = true;
flagFftshift = 'off';

if isempty(s.zeropad)
  zeropad = length(s.time);
else
  zeropad = s.zeropad;
end

if nargin==2
  zeropad = varargin{1};
end
if nargin>=3
  while length(varargin)>=2
    arg = varargin{1};
    val = varargin{2};
    switch lower(arg)
      case {'flagspectrometer','spectrometer'}
        flagSpectrometer = val;
      case {'flagfftshift','fftshift'}
        flagFftshift = val;
      case {'zeropad'}
        zeropad = val;
      otherwise
        error(['unknown option ',arg])
    end
    varargin = varargin(3:end);
  end
end

if zeropad 
  s.zeropad = zeropad;
end
if s.zeropad ~= length(s.time)
  s.comment = [s.comment,' zeropad ',num2str(zeropad)];
end
if isempty(s.undersampling)
  error('in freq2d: undersampling is empty')
end
if s.undersampling<0
  error('in freq2d: undersampling is < 0!')
end
  
[w,p]=fftFreqAxis(s.time,...
  'fftshift',flagFftshift,...
  'freq_units',s.freq_units,...
  'time_units',s.time_units,...
  'undersampling',s.undersampling,...
  'zeropad',zeropad);

s.w1 = w;
if flagSpectrometer==false
  s.w3 = w;
end

s.resolution = p.resolution;
s.centerfreq = p.centerfreq;

if flagSpectrometer
  if ~isempty(s.spec_calib)
    calib = s.spec_calib;
  else
    calib = 0;
  end
  s.w3 = s.freq - calib ;
end
