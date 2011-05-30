function s = absorptive2d(s,varargin)

zeropad = length(s.time);
phase = 0;

fft_type = 'sgrsfft';
fft_type_list = {'fft','sgrsfft'};

apodization = 'none';
apodization_list = {'none','triangular','gaussian'};

flag_plot = true;
n_contours = 20;
range = [1000 3000];

if isa(varargin{1},'char')
  input_arguments_version = 2;
else
  input_arguments_version = 1;
end

switch input_arguments_version
  case 1
    if nargin >= 2
      phase = varargin{1}(1);
      s.phase = phase;
    end
    if nargin >=3
      if ~isempty(varargin{2})
        zeropad = varargin{2};
      end
    end
    if nargin >=4
      if ~isempty(varargin{3})
        range = varargin{3};
      end
    end
    if nargin >=5
      fft_type = varargin{4};
    end
    if nargin >= 6
      if ~isempty(varargin{5})
        apodization = varargin{5};
      end
    end
    
  case 2
    while length(varargin)>=2
      arg = varargin{1};
      val = varargin{2};
      switch lower(arg)
        case 'n_contours'
          n_contours = val;
        case 'phase'
          phase = val(1); %take only the first element if it is an array
        case 'zeropad'
          zeropad = val;
        case 'range'
          range = val;
        case 'fft_type'
          fft_type = val;
        case 'apodization'
          apodization = val;
        %case {'pumpprobe_style','pumpprobe'}
        %  flag_pumpprobe = val;
        case 'plot'
          flag_plot = val;
        otherwise
          error(['rbAbsorptive: unknown option ',arg])
      end
      varargin = varargin(3:end);
    end
end

% check if it is a experimental or simulated spectrum
n_freq = length(s.freq);
if n_freq == 0
  flag_spectrometer = false;  
else
  flag_spectrometer = true;
end
n_time = length(s.time);

% error checking
if mod(n_contours,2)
  warning('myAbsorptive2d: Odd number of contour lines may produce unexpected results!')
end
if ~any(strcmpi(fft_type,fft_type_list))
  error(['fft type ',fft_type,' not known in absorptive2d.m']);
end
if ~any(strcmpi(apodization,apodization_list))
  error(['apodization type ',apodization,' not known in absorptive3d.m']);
end

if flag_spectrometer
  % measured spectrum
  R1 = zeros(zeropad,n_freq);
  R2 = zeros(zeropad,n_freq);

  switch apodization
    case 'none'
      window_fxn = ones(n_time, 1);
    case 'triangular'
      window_fxn = linspace(1,0,n_time)';
    case 'gaussian'
      window_fxn = exp(-(linspace(0,3,n_time)).^2)';
  end
  
  size(s.R1(:,1))
  size(window_fxn)
  
  for i = 1:n_freq
    switch fft_type
      case 'fft'
        R1(:,i) = fft(s.R1(:,i)' .* window_fxn, zeropad);
        R2(:,i) = fft(s.R2(:,i)' .* window_fxn, zeropad);
      case 'sgrsfft'
        R1(:,i) = sgrsfft(s.R1(:,i) .* window_fxn, zeropad);
        R2(:,i) = sgrsfft(s.R2(:,i) .* window_fxn, zeropad);
    end
  end  

  R1 = R1.*exp((phase)*sqrt(-1));
  R2 = R2.*exp(-(phase)*sqrt(-1));
  
  R1 = flipdim(R1,2);
  
  s.R  = real(R1+R2);
  
  ind1 = find(s.w1>=range(1) & s.w1<=range(2));
  ind3 = find(s.w3>=range(1) & s.w3<=range(2));

  my2dPlot(s.w3(ind3),s.w1(ind1),s.R(ind1,ind3), 'n_contours', 20);%, 'zlimit', 35);

  
  
%else
  
  
end


















