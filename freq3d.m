function s = freq3d(s,varargin)
%freq3d set up the frequency axis with undersampling etc
%function s = freq3d(s,[zeropad])
%
%function s = freq3d(s,'arg',val)
%arg can be 
%spectrometer = {|true|,false}
%fftshift = {|'off'|,'on'}
%zeropad = val
%direction = {'fft' 'sgrsfft' 'ifft' 'sgrsifft'}
% fft_type equivalent to 'direction'

flagSpectrometer = true;
flagFftshift = 'off';
zeropad = length(s.time);
fft_type = 'fft';

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
      case {'fft_type'}
        fft_type = val;
      otherwise
        error(['unknown option ',arg])
    end
    varargin = varargin(3:end);
  end
end

if zeropad
  s.zeropad = zeropad;
  s.comment = [s.comment,' zeropad ',num2str(zeropad)];
end
if isempty(s.undersampling)
  error('in freq3d: undersampling is empty')
end
if s.undersampling<0
  error('in freq3d: undersampling is < 0')
end
  
[w,p]=fftFreqAxis(s.time,...
  'fftshift',flagFftshift,...
  'freq_units',s.freq_units,...
  'time_units',s.time_units,...
  'undersampling',s.undersampling,...
  'fft_type',fft_type,...
  'zeropad',zeropad);

s.w1 = w;
s.w3 = w;
s.resolution = p.resolution;
s.centerfreq = p.centerfreq;
if flagSpectrometer==false
  s.w5 = w;
end

if flagSpectrometer
  if isempty(s.spec_calib)
    calib = 0;
  else
    calib = s.spec_calib;
  end
  s.w5 = s.freq-calib;
end
