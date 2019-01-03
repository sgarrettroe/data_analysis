function s = absorptive2d(s,varargin)
%calculate the absorptive spectrum from 2d data
%
%  - phase: use a phase, in radian
%
% - zeropad (number): The zeropadded length. Should be equal to twice the 
%   number of time points for the optimum amount of information in the real 
%   spectrum (the default).
%
% - range ([lim_l lim_u]): Plots over the frequency window of interest 
%   given by the lower and upper limits
%
% - fft_type (name): Type can be 'fft', the normal fft, or 'sgrsfft' which 
%   scales the first data point by 0.5. Default is sgrsfft. 
%
% - apodization (name): Can be 'none', 'triangular', 'gaussian', 'rbOnes',
%   'rbGauss' or 'test'. Others can be implemented by adding the methods to 
%   the apodization_list and then changing the window_fxn. Default is none.
%
% - apod_numbers ([a b]): For the rbOnes, rbGauss and test functions, this
%   is the input needed. For rbOnes and rbGauss it determines the length
%   where the window is set to 0 (a) and how long the Gaussian increase to
%   1 takes (b). For test it is the factor in the exponential (a) and the
%   factor in the gaussian (b). 
%
% - pumpprobe (BOOL): The default behavior is to plot 'pump-probe' style. 
%   Default is true. 
%
% - plot (BOOL): plot the rephasing and non-rephasing spectrum and the 
%   apodization function. Default is true.
%


%default values
flag_debug = false;
n_contours = 20;
zeropad = length(s.time); %means no zeropadding
phase = 0;
range = [2300 2700];
fft_type = 'sgrsfft';
fft_type_list = {'fft','sgrsfft'};
apodization = 'none';
apodization_list = {'none','triangular','gaussian', 'rbOnes', 'rbGauss', 'test'};
apod_numbers = [-0.5 3];
apod_pixel = 16;
flag_pumpprobe = true;
flag_plot = true;
flag_fftshift = 'off';
zeropad = 2*length(s.time);


% for the special window function
apod_numbers = [10 10]; 

%determine which version of the input arguments are being passed based on
%if the first value is a property string or a phase
if isa(varargin{1},'char')
  input_arguments_version = 2;
else
  input_arguments_version = 1;
end

switch input_arguments_version
  case 1
    if nargin >= 2
      phase = varargin{1}(1);
      %phase2 = varargin{1}(2);
      s.phase = phase;
      %s.phase2 = phase2;
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
          if mod(n_contours,2)
            warning('my2dPlot: Odd number of contour lines may produce unexpected results!')
          end
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
        case 'apod_numbers'
          apod_numbers = val;
        case 'apod_pixel'
          apod_pixel = val;
        case {'pumpprobe_style','pumpprobe'}
          flag_pumpprobe = val;
        case 'plot'
          flag_plot = val;
        case 'debug'
          flag_debug = val;
        otherwise
          error(['my2dPlot: unknown option ',arg])
      end
      varargin = varargin(3:end);
    end
end

n_freq = length(s.freq);
if n_freq == 0
  flag_spectrometer = false;
  flag_remove_DC=false;
  flag_plot=false;
  if flag_debug; disp(['absorptive2d: flag_spectrometer = false']); end
else
  flag_spectrometer = true;
  flag_remove_DC=true;
  if flag_debug; disp(['absorptive2d: flag_spectrometer = true']); end
end
n_time = length(s.time);

if flag_debug; disp(['absorptive2d: n_time:' int2str(n_time) ', n_freq:' int2str(n_freq)]); end

%error checking of inputs here?
if ~any(strcmpi(fft_type,fft_type_list)), error(['fft type ',fft_type,' not known in absorptive2d.m']);end
if ~any(strcmpi(apodization,apodization_list)), error(['apodization type ',apodization,' not known in absorptive3d.m']);end
s.comment = [s.comment,' // fft-type ',fft_type,' apodization ', apodization];

if flag_spectrometer
  %begin calculation
  R1 = zeros(zeropad,n_freq);
  R2 = zeros(zeropad,n_freq);
  
  switch apodization
    case 'none'
      window_fxn = ones(1, n_time);
      if flag_debug == true
        figure(1001),clf;
        plot(s.R1(:, apod_pixel) .* window_fxn');
      end
    case 'triangular'
      window_fxn = linspace(1,0,n_time);
      if flag_debug == true
        figure(1000);
        plot(window_fxn);         
        figure(1001),clf;
        plot(s.R1(:, apod_pixel) .* window_fxn');
      end
    case 'gaussian'
      window_fxn = exp(-(linspace(0,3,n_time)).^2);
      if flag_debug == true
        figure(1000);
        plot(window_fxn);
        figure(1001),clf;
        plot(s.R1(:, apod_pixel) .* window_fxn');
      end
    case 'rbOnes'
      % Gaussian
      number_a = apod_numbers(1);
      number_b = apod_numbers(2);
      
      a = zeros(1, number_a);
      %b = 1/sqrt(pi * std) * exp(-linspace(-1, 0, number_b).^2 / std) - 0.039;
      b = exp(-(5/1)*linspace(-1, 0, number_b).^2);
      c = ones(1, n_time - number_a - number_b);
      
      window_fxn = cat(2, a, b, c);

      if flag_plot == true
        figure(1000);
        plot(window_fxn);
        figure(1001);
        plot(s.R1(:, apod_pixel) .* window_fxn');
      end

    case 'rbGauss'
      % Gaussian
      number_a = apod_numbers(1);
      number_b = apod_numbers(2);
      
      a = zeros(1, number_a);
      b = exp(-(5/1)*linspace(-1, 0, number_b).^2);
      c = exp(-(linspace(0, 3, n_time - number_a - number_b)).^2);
      
      window_fxn = cat(2, a, b, c);

      if flag_plot == true
        figure(1000);
        plot(window_fxn);
        figure(1001);
        plot(s.R1(:, apod_pixel) .* window_fxn');
      end
      
    case 'test'
      number_a = apod_numbers(1);
      number_b = apod_numbers(2);
      
      %a = zeros(1, number_a);
      b = exp(- linspace(0, number_a, n_time)) - 1;
      c = exp(-(linspace(0, number_b, n_time)).^2);
      
      window_fxn = b .* c; % cat(2, b, c);
      
      window_fxn = window_fxn / max(window_fxn);
      
      if flag_plot == true
        figure(1000),clf;
        hold on;
        plot(window_fxn);
        plot(b);
        plot(c);
        hold off;
        figure(1001),clf;
        plot(s.R1(:, apod_pixel) .* window_fxn');
      end      
  end
  % end switch apodization
  
  for i = 1:n_freq
    switch fft_type
      case 'fft'
        R1(:,i) = fft(s.R1(:,i)',zeropad);
        R2(:,i) = fft(s.R2(:,i)',zeropad);
      case 'sgrsfft'
        R1(:,i) = sgrsfft(s.R1(:,i) .* window_fxn', zeropad);
        R2(:,i) = sgrsfft(s.R2(:,i) .* window_fxn', zeropad);
    end
  end
  
  R1 = R1.*exp((phase)*sqrt(-1));
  R2 = R2.*exp(-(phase)*sqrt(-1));
  
  if flag_remove_DC
    R1(1,:)=0;
    R2(1,:)=0;
  end
  
  %do the flips
  flag_flips = false;
  if flag_flips
    if mod(length(s.R2),2)==0
      %if an even number of points be careful about the zero frequency
      R1 = fliplr(circshift(R1,[0 -1]));
      %R1 = flipdim(circshift(R1,[0 -1 0]),2);
    else
      %this needs to be double checked...
      %R1 = fliplr(circshift(R1,[0 -1]));
      R1 = flipdim(R1,2);
    end
  end
  s.R  = real(R1+R2);
  
  %redo frequency axis in case we zeropadded
  %s = freq2d(s,zeropad);
  s = freq2d(s,'zeropad',zeropad,...
    'spectrometer',flag_spectrometer,...
    'fftshift',flag_fftshift);
  
  
  
  
  map = myMapRGB2(n_contours+1);
  ind = find(s.w1>range(1) & s.w1<range(2));
  
  if flag_plot
    figure(100)
    subplot(1,2,1)
    if flag_pumpprobe
      x = s.freq;
      y = s.w1(ind);
      z = real(R1(ind,:));
      x_label = '\omega_{probe} / 2\pic';
      y_label = '\omega_{pump} / 2\pic'; 
    else
      x = s.w1(ind);
      y = s.freq;
      z = real(R1(ind,:)');
      x_label = '\omega_1 / 2\pic';
      y_label = '\omega_3 / 2\pic'; 
    end
    
    contourf(x,y,z,n_contours)
    axis square
    myCaxis2(z,n_contours);
    colormap(map)
    xlabel(x_label)
    ylabel(y_label)
    
    subplot(1,2,2)
    if flag_pumpprobe
      x = s.freq;
      y = s.w1(ind);
      z = real(R2(ind,:));
      x_label = '\omega_{probe} / 2\pic';
      y_label = '\omega_{pump} / 2\pic'; 
    else
      x = s.w1(ind);
      y = s.freq;
      z = real(R1(ind,:)');
      x_label = '\omega_1 / 2\pic';
      y_label = '\omega_3 / 2\pic'; 
    end
    contourf(x,y,z,n_contours)
    axis square
    myCaxis2(z,n_contours);
    colormap(map)
    xlabel(x_label)
    ylabel(y_label)
    
    figure(101),clf
    if flag_pumpprobe
      x = s.w3;
      y = s.w1(ind);
      z = s.R(ind,:);
      x_label = '\omega_{probe} / 2\pic';
      y_label = '\omega_{pump} / 2\pic'; 
    else
      x = s.w1(ind);
      y = s.w3;
      z = s.R(ind,:)';
      x_label = '\omega_1 / 2\pic';
      y_label = '\omega_3 / 2\pic'; 
    end
     a=my2dPlot(x,y,z,'n_contours',n_contours,'pumpprobe',flag_pumpprobe);
  end %if flag_plot
  
else
  %if time domain experiment
  %begin calculation
  R1 = zeros(zeropad,zeropad);
  R2 = zeros(zeropad,zeropad);
    
  switch apodization
    case 'none'
      window_fxn = ones(n_time,n_time);
    case 'triangular'
      temp = linspace(1,0,n_time);
      [X,Y] = meshgrid(temp,temp);
      window_fxn = X.*Y;
    case 'gaussian'
      temp = linspace(0,3,n_time);
      [X,Y] = meshgrid(temp,temp);
      window_fxn = exp(-(X+Y).^2);
  end
 
  switch fft_type
    case 'fft'
      R1 = fftn(s.R1.*window_fxn,[zeropad,zeropad]);
      R2 = fftn(s.R2.*window_fxn,[zeropad,zeropad]);
    case 'sgrsfft'
      R1 = sgrsfft2(s.R1.*window_fxn,zeropad);
      R2 = sgrsfft2(s.R2.*window_fxn,zeropad);
  end
  
  R1 = R1.*exp((phase)*sqrt(-1));
  R2 = R2.*exp(-(phase)*sqrt(-1));
  
  if flag_remove_DC
    R1(1,:)=0;
    R2(1,:)=0;
  end
  
  %do the flips
  flag_flips = false;
  if flag_flips
    if mod(length(s.R2),2)==0
      %if an even number of points be careful about the zero frequency
      R1 = fliplr(circshift(R1,[0 -1]));
      %R1 = flipdim(circshift(R1,[0 -1 0]),2);
    else
      %this needs to be double checked...
      %R1 = fliplr(circshift(R1,[0 -1]));
      R1 = flipdim(R1,2);
    end
  end
  s.R  = real(R1+R2);
  
  %redo frequency axis if we zeropadded
  if zeropad~=length(s.time)
    s = freq2d(s,zeropad);
  end
   
  %if plotting
  if flag_plot
     figure(101),clf
     
    if flag_pumpprobe
      x = s.w3;
      y = s.w1;
      z = s.R'; %not sure about this
      x_label = '\omega_{probe} / 2\pic';
      y_label = '\omega_{pump} / 2\pic'; 
    else
      x = s.w1;
      y = s.w3;
      z = s.R; %not sure about this 
      x_label = '\omega_1 / 2\pic';
      y_label = '\omega_3 / 2\pic'; 
    end
    a = my2dPlot(x,y,z,'n_contours',n_contours,'pumpprobe',flag_pumpprobe);
  end %if flag_plot

end
