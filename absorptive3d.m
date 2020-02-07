function s = absorptive3d(s,varargin)
%calculate the absorptive spectrum from 3d data
%syntax: s = absorptive3d(s,[phase,ipix,zeropad,fft_type,apodization,range,fftshift])
n_freq = length(s.freq);
if n_freq == 0
  flagSpectrometer = false;
  flagRemoveDC = false;
else
  flagSpectrometer = true;
  flagRemoveDC = true;
end
n_time = length(s.time);

flagPlot = false;
flag3dPlot = true;
flagFftshift = 'off';

%default values
zeropad = n_time; %means no zeropadding
phase = 0;
ipix = 16;
fft_type= 'fft';
fft_type_list = {'fft','sgrsfft','ifft','sgrsifft'};
apodization = 'none';
apodization_list = {'none','triangular','gaussian'};
range = [2300 2700];
if nargin >= 2
  phase = varargin{1};
  s.phase = phase;
  switch length(phase)
    case 1
      flag_phase_differences = false;
    case 4
      flag_phase_differences = true;
    otherwise
      error(['absorptive.m phase array should have 1 or 4 elements, it has ', num2str(length(phase))]);
  end
end
if nargin >=3
  ipix = varargin{2};
end
if nargin >=4
  if ~isempty(varargin{3})
    zeropad = varargin{3};
  end
end
if nargin >= 5
  if ~isempty(varargin{4})
    fft_type = varargin{4};
  end
end
if nargin >= 6
  if ~isempty(varargin{5})
    apodization = varargin{5};
  end
end
if nargin >=7
  if ~isempty(varargin{6})
    range = varargin{6};
  end
end
if nargin >=8
  if ~isempty(varargin{7})
    flagFftshift = varargin{7};
  end
end

%error checking of inputs here?
fft_type = lower(fft_type);
if ~any(strcmpi(fft_type,fft_type_list)), error(['fft type ',fft_type,' not known in absorptive3d.m']);end
if ~any(strcmpi(apodization,apodization_list)), error(['apodization type ',apodization,' not known in absorptive3d.m']);end
s.comment = [s.comment,' fft_type ',fft_type,' apodization ',apodization];

%------------------------------------------------------------------------
%
%     Do the FFTs
%
%------------------------------------------------------------------------
if flagSpectrometer
  R1 = zeros(zeropad,zeropad,n_freq);
  R2 = zeros(zeropad,zeropad,n_freq);
  R3 = zeros(zeropad,zeropad,n_freq);
  R4 = zeros(zeropad,zeropad,n_freq);

  switch apodization
    case 'none'
      window_fxn = ones(n_time,n_time);
    case 'triangular'
      [X,Y] = meshgrid(linspace(1,0,n_time));
      window_fxn = X.*Y;
    case 'gaussian'
      [X,Y] = meshgrid(linspace(0,3,n_time));
      window_fxn = exp(-(X+Y).^2);
  end
  
  switch fft_type
    case 'fft'
      for i = 1:n_freq
        R1(:,:,i) = fft2(squeeze(s.R1(:,:,i)).*window_fxn,zeropad,zeropad);
        R2(:,:,i) = fft2(squeeze(s.R2(:,:,i)).*window_fxn,zeropad,zeropad);
        R3(:,:,i) = fft2(squeeze(s.R3(:,:,i)).*window_fxn,zeropad,zeropad);
        R4(:,:,i) = fft2(squeeze(s.R4(:,:,i)).*window_fxn,zeropad,zeropad);
      end
    case 'sgrsfft'
      for i = 1:n_freq
        R1(:,:,i) = sgrsfft2(squeeze(s.R1(:,:,i)).*window_fxn,zeropad,zeropad);
        R2(:,:,i) = sgrsfft2(squeeze(s.R2(:,:,i)).*window_fxn,zeropad,zeropad);
        R3(:,:,i) = sgrsfft2(squeeze(s.R3(:,:,i)).*window_fxn,zeropad,zeropad);
        R4(:,:,i) = sgrsfft2(squeeze(s.R4(:,:,i)).*window_fxn,zeropad,zeropad);
      end
    case 'ifft'
      for i = 1:n_freq
        R1(:,:,i) = ifft2(squeeze(s.R1(:,:,i)).*window_fxn,zeropad,zeropad);
        R2(:,:,i) = ifft2(squeeze(s.R2(:,:,i)).*window_fxn,zeropad,zeropad);
        R3(:,:,i) = ifft2(squeeze(s.R3(:,:,i)).*window_fxn,zeropad,zeropad);
        R4(:,:,i) = ifft2(squeeze(s.R4(:,:,i)).*window_fxn,zeropad,zeropad);
      end
    case 'sgrsifft'
      for i = 1:n_freq
        R1(:,:,i) = sgrsifft2(squeeze(s.R1(:,:,i)).*window_fxn,zeropad,zeropad);
        R2(:,:,i) = sgrsifft2(squeeze(s.R2(:,:,i)).*window_fxn,zeropad,zeropad);
        R3(:,:,i) = sgrsifft2(squeeze(s.R3(:,:,i)).*window_fxn,zeropad,zeropad);
        R4(:,:,i) = sgrsifft2(squeeze(s.R4(:,:,i)).*window_fxn,zeropad,zeropad);
      end
  end

else %if time domain experiment / simulation
  R1 = zeros(zeropad,zeropad,zeropad);
  R2 = zeros(zeropad,zeropad,zeropad);
  R3 = zeros(zeropad,zeropad,zeropad);
  R4 = zeros(zeropad,zeropad,zeropad);
  
  switch apodization
    case 'none'
      window_fxn = ones(n_time,n_time,n_time);
    case 'triangular'
      temp = linspace(1,0,n_time);
      [X,Y,Z] = meshgrid(temp,temp,temp);
      window_fxn = X.*Y.*Z;
    case 'gaussian'
      temp = linspace(0,3,n_time);
      [X,Y,Z] = meshgrid(temp,temp,temp);
      window_fxn = exp(-(X+Y+Z).^2);
  end
 
  switch fft_type
    case 'fft'
      R1 = fftn(s.R1.*window_fxn,[zeropad,zeropad,zeropad]);
      R2 = fftn(s.R2.*window_fxn,[zeropad,zeropad,zeropad]);
      R3 = fftn(s.R3.*window_fxn,[zeropad,zeropad,zeropad]);
      R4 = fftn(s.R4.*window_fxn,[zeropad,zeropad,zeropad]);
    case 'sgrsfft'
      R1 = sgrsfft3(s.R1.*window_fxn,zeropad);
      R2 = sgrsfft3(s.R2.*window_fxn,zeropad);
      R3 = sgrsfft3(s.R3.*window_fxn,zeropad);
      R4 = sgrsfft3(s.R4.*window_fxn,zeropad);
    case 'ifft'
      R1 = ifftn(s.R1.*window_fxn,[zeropad,zeropad,zeropad]);
      R2 = ifftn(s.R2.*window_fxn,[zeropad,zeropad,zeropad]);
      R3 = ifftn(s.R3.*window_fxn,[zeropad,zeropad,zeropad]);
      R4 = ifftn(s.R4.*window_fxn,[zeropad,zeropad,zeropad]);
    case 'sgrsifft'
      R1 = sgrsifft3(s.R1.*window_fxn,zeropad);
      R2 = sgrsifft3(s.R2.*window_fxn,zeropad);
      R3 = sgrsifft3(s.R3.*window_fxn,zeropad);
      R4 = sgrsifft3(s.R4.*window_fxn,zeropad);
  end

end

%redo frequency axis if we zeropadded
s = freq3d(s,'zeropad',zeropad,...
  'spectrometer',flagSpectrometer,...
  'fftshift',flagFftshift,...
  'fft_type',fft_type);
  
%------------------------------------------------------------------------
%
%     Phase and add spectra
%
%------------------------------------------------------------------------
if flag_phase_differences
  R1 = R1.*exp(1i*phase(1));
  R2 = R2.*exp(1i*phase(2));
  R3 = R3.*exp(1i*phase(3));
  R4 = R4.*exp(1i*phase(4));
else
  R1 = R1.*exp(1i*phase);
  R2 = R2.*exp(1i*phase);
  R3 = R3.*exp(1i*phase);
  R4 = R4.*exp(1i*phase);
end

if flagSpectrometer
  if flagRemoveDC
    R1(1,1,:)=0;
    R2(1,1,:)=0;
    R3(1,1,:)=0;
    R4(1,1,:)=0;
  end
  
  R2 = flipdim(circshift(R2,[0 -1 0]),2);
  R3 = flipdim(circshift(R3,[-1 0 0]),1);
  R4 = flipdim(flipdim(circshift(R4,[-1 -1 0]),1),2);
  s.R  = real(R1+R2+R3+R4);
  
  
  map = myMapRGB(65);
  ind = find(s.w1>range(1) & s.w1<range(end));
  if isempty(ind)
    warning('sgr:emptyRange', ...
      'Found no data in the given frequency range')
    ind = 1:length(s.w1);
  end
  
  if flagPlot
    figure(100),clf
    subplot(2,2,1)
    contourf(s.w1(ind),s.w3(ind),real(R1(ind,ind,ipix)))
    axis square
    myCaxis(real(R1(ind,ind,ipix)));
    colormap(map)
    xlabel('\omega_1 / 2\pic')
    ylabel('\omega_3 / 2\pic')
    
    subplot(2,2,2)
    contourf(s.w1(ind),s.w3(ind),real(R2(ind,ind,ipix)))
    axis square
    myCaxis(real(R2(ind,ind,ipix)));
    colormap(map)
    xlabel('\omega_1 / 2\pic')
    ylabel('\omega_3 / 2\pic')
    
    subplot(2,2,3)
    contourf(s.w1(ind),s.w3(ind),real(R3(ind,ind,ipix)))
    axis square
    myCaxis(real(R3(ind,ind,ipix)));
    colormap(map)
    xlabel('\omega_1 / 2\pic')
    ylabel('\omega_3 / 2\pic')
    
    subplot(2,2,4)
    contourf(s.w1(ind),s.w3(ind),real(R4(ind,ind,ipix)))
    axis square
    myCaxis(real(R4(ind,ind,ipix)));
    colormap(map)
    xlabel('\omega_1 / 2\pic')
    ylabel('\omega_3 / 2\pic')
    
    figure(101),clf
    contourf(s.w1(ind),s.w3(ind),s.R(ind,ind,ipix),10);
    axis square
    colormap(map)
    myCaxis(s.R(ind,ind,ipix));
    xlabel('\omega_1 / 2\pic')
    ylabel('\omega_3 / 2\pic')
    
    figure(103),clf
    subplot(2,2,1)
    my3dPlot3(s.w1(ind),s.w3(ind),s.w5,real(R1(ind,ind,:)),'clf','off')
    subplot(2,2,2)
    my3dPlot3(s.w1(ind),s.w3(ind),s.w5,real(R2(ind,ind,:)),'clf','off','zlabel','off')
    subplot(2,2,3)
    my3dPlot3(s.w1(ind),s.w3(ind),s.w5,real(R3(ind,ind,:)),'clf','off')
    subplot(2,2,4)
    my3dPlot3(s.w1(ind),s.w3(ind),s.w5,real(R4(ind,ind,:)),'clf','off','zlabel','off')
    
  end %if flagPlot
  
  if flag3dPlot
    figure(102),clf
    my3dPlot(s.w1(ind),s.w3(ind),s.w5,s.R(ind,ind,:))
  end %if flag3dPlot
  
else %if all time domain
  
  %do the flips
  R2 = flipdim(circshift(R2,[0 -1 0]),2);
  R3 = flipdim(circshift(R3,[-1 0 0]),1);
  R4 = flipdim(flipdim(circshift(R4,[-1 -1 0]),1),2);
  s.R = fftshift(real(R1+R2+R3+R4));
  
end
