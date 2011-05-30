function s = absorptive2d(s,varargin)
%calculate the absorptive spectrum from 2d data
%
% s = absorptive2d(s,'Property',value,...)
%
% s = absorptive2d(s,'phase',val)
%     uses a phase of val
%
% s = absorptive2d(s,'zeropad',val)
%     The zeropadded length. Should be equal to twice the number of time
%     points for the optimum amount of information in the real spectrum
%     (the default).
%
% s = absorptive2d(s,'range',[lim_l lim_u])
%     Plots over the frequency window of interest given by the lower and 
%     upper limits
%
% s = absorptive2d(s,'fft_type','type')
%     Type can be 'fft', the normal fft, or 'sgrsfft' which scales the 
%     first data point by 0.5 
%
% s = absorptive2d(s,'apodization','type')
%     Can be triangular or gaussian. Others can be implemented by adding
%     the methods to the apodization_list and then changing the window_fxn
%
% s = absorptive2d(s,'pumpprobe',true) 
%     The default behavior is to plot 'pump-probe' style 
%
% s = absorptive2d(s,'pumpprobe',false) 
%     The plots are (x,y) = (omega_1, omega_3)  style 
%
% s = absorptive2d(s,'plot',true) 
% s = absorptive2d(s,'plot',false) 
%     Turn the plots on or off
%
% 

%default values
n_contours = 20;
zeropad = length(s.time); %means no zeropadding
phase = 0;
range = [2300 2700];
fft_type = 'sgrsfft';
fft_type_list = {'fft','sgrsfft'};
apodization = 'none';
apodization_list = {'none','triangular','gaussian'};
flag_pumpprobe = true;
flag_plot=true;
flag_fftshift = 'off';
zeropad = 2*length(s.time);

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
        case {'pumpprobe_style','pumpprobe'}
          flag_pumpprobe = val;
        case 'plot'
          flag_plot = val;
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
else
  flag_spectrometer = true;
  flag_remove_DC=true;
end
n_time = length(s.time);



%error checking of inputs here?
if ~any(strcmpi(fft_type,fft_type_list)), error(['fft type ',fft_type,' not known in absorptive2d.m']);end
if ~any(strcmpi(apodization,apodization_list)), error(['apodization type ',apodization,' not known in absorptive3d.m']);end
s.comment = [s.comment,' fft_type ',fft_type,' apodization ',apodization];

if flag_spectrometer
  %begin calculation
  R1 = zeros(zeropad,n_freq);
  R2 = zeros(zeropad,n_freq);
  for i = 1:n_freq
    switch fft_type
      case 'fft'
        R1(:,i) = fft(s.R1(:,i)',zeropad);
        R2(:,i) = fft(s.R2(:,i)',zeropad);
      case 'sgrsfft'
        R1(:,i) = sgrsfft(s.R1(:,i),zeropad);
        R2(:,i) = sgrsfft(s.R2(:,i),zeropad);
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
  figure(100),
  
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
