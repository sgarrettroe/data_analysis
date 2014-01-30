function [out,extra] = analyticalResponseFunctionsFun(p,w1_in,w3_in,options)
%This is an example of how to use analyticalResponseFunctions.m
%all time parameters are in units of ps

%Here are the things you need to define:
flag_print = 0; %1 => figures or 0 => no figures
flag_plot = 0;
order = 3; %order of spectroscopy to calculate. 3 = 2DIR, 5 = 3DIR

%----------------------
%
%  system properties 
%
%----------------------
%specify the form and parameters of the frequency fluctuation correlation
%functions. Uncomment a block to use that form and its parameters

%%for a single overdamped motion 
%damping = 'overdamped'; %i.e. exponential decay
%Delta_cm = 10; %linewidth (sigma) in wavenumbers
%tau = 2.; % correlation time in ps

%%for critical damping (slightly more oscillation in the correlation function)
%damping = 'critical';
%Delta_cm = 10; %linewidth (sigma) in wavenumbers
%tau = 2.; % correlation time in ps

%%for water one can use a fit to results from simulation (all parameters are
%%hard coded basically)
%damping = 'hynesform';

%for multiexponential decay (fast and slow motion). 
%damping = 'multiexp';

damping = options.damping;
% Delta1_cm = 8;%linewidth (sigma) in wavenumbers of one motion
% Delta2_cm = 8;%linewidth (sigma) in wavenumbers of the other motion
% tau1 = 0.2; %first timescale (ps)
% tau2 = 15; %second timescale (ps)

%population times to calculate
%t2_array = 0.2; %ps (can be a single time)
%t2_array = [0.2 10 50]; %can be a vector of times
t2_array = options.t2_array; %can be a vector of times

%for 5th order also define t4 times (not used for 2DIR)
t4_array = t2_array;
n_t2_array = length(t2_array);

% THIS NEEDS FIXING!!!

%for an exponential c2(t)
%dt = 2*0.0021108; %ps
%dt = 0.02;
%n_t = 256; %number of time steps
dt = options.dt;
n_t = options.n_t; %number of time steps
w_0_cm = options.w_0_cm;% %center frequency
if isfield(options,'order')
    order = options.order;
end

flag_rotating_frame = true;

range = [1960 2160];
labels = [-400:200:400]; %not used in 2d plots
%range = [-200 200];
%labels = [-150:150:150];

two_level_system = false;
%two_level_system = true;


%details of fft
%fft_type = 'fft';
%fft_type = 'petersfft';
fft_type = 'sgrsfft';
n_interp = 1; %number of points of linear interpolation to use
n_zp = 2*n_t; %total length after zeropadding 
            %(make n_zp=n_t for no zero padding)
n_under = 0;%2;

apodization = 'none';
%apodization = 'triangular';
%apodization = 'gaussian';

%type of projection to calculate
%projection_type = 'window';
projection_type = 'all';

% play with phase shifts
phi = 0; %phase shift in deg
noise = 0.0;

% simulate laser bandwidth
simulate_bandwidth = false;
%simulate_bandwidth = true;
bandwidth = 180; %cm-1 fwhm
bandwidth_axes = 2; %2 for scaled by LO spectrum, 3 for without

%
% Body of the calculation
%
extra = [];

%-------------------------------------------------------------
%
% start calculation
%
%-------------------------------------------------------------
c = 2.9979e10;
wavenumbersToInvPs=c*1e-12;
invPsToWavenumbers=1/wavenumbersToInvPs;
c_cmfs = c*1e-15;
t=0:dt:(n_t-1)*dt;

w = fftFreqAxis(t,...
  'time_units','ps',...
  'freq','wavenumbers',...
  'shift','on',...
  'zeropad',n_zp,...
  'undersampling',n_under);
if flag_rotating_frame
    w = w + w_0_cm;
end

dw = w(2)-w(1);
w_0 = w_0_cm*2*pi*wavenumbersToInvPs; %convert to radians

if ~exist('orientational_response','var')
  %if it is not defined then assume you don't care
  orientational_response = false;
end
% if orientational_response
%   disp('Using orentational response functions');
% else
%   disp('No orentational response functions');
% end
  
c2 = [];
g = [];
switch damping,
    case 'voigt'
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    T2 = p(2); %first timescale (ps)
    anh_cm = p(3);
        
    Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;

    c2 = @(t) (t==0)/T2 + Delta1^2;
    g = @(t) t./T2 + Delta1^2/2.*t.^2;
        
    case {'overdamped', '1exp'}
    %overdamped exp(-t/tau)
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    tau1 = p(2); %first timescale (ps)
    anh_cm = p(3);

    Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
    Lambda1 = 1/tau1;
    anh = anh_cm*wavenumbersToInvPs*2*pi;

    c2 = @(t) Delta1^2.*exp(-Lambda1.*t);
    g = @(t) Delta1^2/Lambda1^2.*(exp(-Lambda1.*t)-1+Lambda1.*t);
  case 'critical'
    %critically damped (1+2t/tau)exp(-2t/tau)
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    tau1 = p(2); %first timescale (ps)
    anh_cm = p(3);

    Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
    Lambda1 = 1/tau1;
    anh = anh_cm*wavenumbersToInvPs*2*pi;

    c2 = @(t) Delta1^2.*(1+2*Lambda1.*t).*exp(-2*Lambda1.*t);
    g = @(t) Delta1^2/4/Lambda1^2.*exp(-2.*Lambda1.*t) ...
      .*(3 + 2*Lambda1*t + exp(2.*Lambda1.*t).*(4*Lambda1.*t - 3));
  case '2expcrit'
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    Delta2_cm = p(2);%linewidth (sigma) in wavenumbers of the other motion
    tau1 = p(3); %first timescale (ps)
    tau2 = p(4); %second timescale (ps)
    anh_cm = p(5);

    Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
    Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
    Lambda1 = 1/tau1;
    Lambda2 = 1/tau2;
    anh = anh_cm*wavenumbersToInvPs*2*pi;
    g = @(t) Delta1^2/4/Lambda1^2 ...
      .*(3.*exp(-2.*Lambda1.*t) + 2*Lambda1*t.*exp(-2.*Lambda1.*t) + (4*Lambda1.*t - 3)) ...
      + Delta2^2/4/Lambda2^2 ...
      .*(3.*exp(-2.*Lambda2.*t) + 2*Lambda2*t.*exp(-2.*Lambda2.*t) + (4*Lambda2.*t - 3));

  case '3expcrit'
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    Delta2_cm = p(2);%linewidth (sigma) in wavenumbers of the other motion
    Delta3_cm = p(3);%linewidth (sigma) in wavenumbers of the other motion
    tau1 = p(4); %first timescale (ps)
    tau2 = p(5); %second timescale (ps)
    tau3 = p(6); %second timescale (ps)
    anh_cm = p(7);

    Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
    Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
    Delta3 = Delta3_cm*wavenumbersToInvPs*2*pi;
    Lambda1 = 1/tau1;
    Lambda2 = 1/tau2;
    Lambda3 = 1/tau3;
    anh = anh_cm*wavenumbersToInvPs*2*pi;
    g = @(t) Delta1^2/4/Lambda1^2 ...
      .*(3.*exp(-2.*Lambda1.*t) + 2*Lambda1*t.*exp(-2.*Lambda1.*t) + (4*Lambda1.*t - 3)) ...
      + Delta2^2/4/Lambda2^2 ...
      .*(3.*exp(-2.*Lambda2.*t) + 2*Lambda2*t.*exp(-2.*Lambda2.*t) + (4*Lambda2.*t - 3)) ...      
      + Delta3^2/4/Lambda3^2 ...
      .*(3.*exp(-2.*Lambda3.*t) + 2*Lambda3*t.*exp(-2.*Lambda3.*t) + (4*Lambda3.*t - 3));

  case '3exp1offcrit'
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    Delta2_cm = p(2);%linewidth (sigma) in wavenumbers of the other motion
    Delta3_cm = p(3);%linewidth (sigma) in wavenumbers of the other motion
    Delta4_cm = p(4);%linewidth (sigma) in wavenumbers of static component
    tau1 = p(5); %first timescale (ps)
    tau2 = p(6); %second timescale (ps)
    tau3 = p(7); %second timescale (ps)
    anh_cm = p(8);

    Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
    Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
    Delta3 = Delta3_cm*wavenumbersToInvPs*2*pi;
    Delta4 = Delta4_cm*wavenumbersToInvPs*2*pi;
    Lambda1 = 1/tau1;
    Lambda2 = 1/tau2;
    Lambda3 = 1/tau3;
    anh = anh_cm*wavenumbersToInvPs*2*pi;
    g = @(t) Delta1^2/4/Lambda1^2 ...
      .*(3.*exp(-2.*Lambda1.*t) + 2*Lambda1*t.*exp(-2.*Lambda1.*t) + (4*Lambda1.*t - 3)) ...
      + Delta2^2/4/Lambda2^2 ...
      .*(3.*exp(-2.*Lambda2.*t) + 2*Lambda2*t.*exp(-2.*Lambda2.*t) + (4*Lambda2.*t - 3)) ...      
      + Delta3^2/4/Lambda3^2 ...
      .*(3.*exp(-2.*Lambda3.*t) + 2*Lambda3*t.*exp(-2.*Lambda3.*t) + (4*Lambda3.*t - 3)) ...
      + Delta4^2.*t.^2/2;

  case '2exp'
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    Delta2_cm = p(2);%linewidth (sigma) in wavenumbers of the other motion
    tau1 = p(3); %first timescale (ps)
    tau2 = p(4); %second timescale (ps)
    anh_cm = p(5);

    Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
    Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
    Lambda1 = 1/tau1;
    Lambda2 = 1/tau2;
    anh = anh_cm*wavenumbersToInvPs*2*pi;
    g = @(t) Delta1^2/Lambda1^2 ...
      .*(exp(-Lambda1.*t) - 1 + Lambda1*t) ...
      + Delta2^2/Lambda2^2 ...
      .*(exp(-Lambda2.*t) - 1 + Lambda2*t);

  case '3exp'
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    Delta2_cm = p(2);%linewidth (sigma) in wavenumbers of the other motion
    Delta3_cm = p(3);%linewidth (sigma) in wavenumbers of the other motion
    tau1 = p(4); %first timescale (ps)
    tau2 = p(5); %second timescale (ps)
    tau3 = p(6); %second timescale (ps)
    anh_cm = p(7);

    Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
    Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
    Delta3 = Delta3_cm*wavenumbersToInvPs*2*pi;
    Lambda1 = 1/tau1;
    Lambda2 = 1/tau2;
    Lambda3 = 1/tau3;
    anh = anh_cm*wavenumbersToInvPs*2*pi;
    g = @(t)Delta1^2/Lambda1^2 ...
      .*(exp(-Lambda1.*t) - 1 + Lambda1*t) ...
      + Delta2^2/Lambda2^2 ...
      .*(exp(-Lambda2.*t) - 1 + Lambda2*t) ...
      + Delta3^2/Lambda3^2 ...
      .*(exp(-Lambda3.*t) - 1 + Lambda3*t);

  case '3exp1off'
    Delta1_cm = p(1);%linewidth (sigma) in wavenumbers of one motion
    Delta2_cm = p(2);%linewidth (sigma) in wavenumbers of the other motion
    Delta3_cm = p(3);%linewidth (sigma) in wavenumbers of the other motion
    Delta4_cm = p(4);%linewidth (sigma) in wavenumbers of static component
    tau1 = p(5); %first timescale (ps)
    tau2 = p(6); %second timescale (ps)
    tau3 = p(7); %second timescale (ps)
    anh_cm = p(8);

    Delta1 = Delta1_cm*wavenumbersToInvPs*2*pi;
    Delta2 = Delta2_cm*wavenumbersToInvPs*2*pi;
    Delta3 = Delta3_cm*wavenumbersToInvPs*2*pi;
    Delta4 = Delta4_cm*wavenumbersToInvPs*2*pi;
    Lambda1 = 1/tau1;
    Lambda2 = 1/tau2;
    Lambda3 = 1/tau3;
    anh = anh_cm*wavenumbersToInvPs*2*pi;
    g = @(t)Delta1^2/Lambda1^2 ...
      .*(exp(-Lambda1.*t) - 1 + Lambda1*t) ...
      + Delta2^2/Lambda2^2 ...
      .*(exp(-Lambda2.*t) - 1 + Lambda2*t) ...
      + Delta3^2/Lambda3^2 ...
      .*(exp(-Lambda3.*t) - 1 + Lambda3*t) ...
      + Delta4^2.*t.^2/2;

  case 'hynesform'
  %fit c2 to hynes type function, double integrate with Mathematica
  %to find g(t)
  a1 = 0.3232;
  k11 = 30.01; %ps-1
  k12 = 17.41; %ps-1
  a2 = 0.3378;
  k2 = 8.270; %ps-1
  a3 = 0.3455;
  k3 = 1.897; %ps-1
  
  %yuck
  g = @(t) Delta^2*exp(-t.*(k12+k2+k3))/(k2^2*k3^2*(k11^2+k12^2)^2) ...
      .*(a1*k2^2*k3^2.*exp((k2+k3).*t).*(cos(k11.*t).*(k12^2-k11^2)-2*k11*k12*sin(k11*t) ...
					 + exp(k12.*t).*(k11^2*(k12.*t+1)+k12^2*(k12.*t-1))) ... 
	 + (k11^2+k12^2)^2.*exp(k12.*t).*(a3*k2^2*exp(k2.*t).*(exp(k3.*t).*(k3.*t-1)+1) ...
					  +a2*k3^2*exp(k3.*t).*(exp(k2.*t).*(k2*t-1)+1)));
  otherwise
    error('damping value is unknown');
end

if order==1
  %S = exp(-g(t)).*cos(w_0.*t);
  S1 = exp(-g(t));
  %S2 = exp(-g(t)).*cos(w_0.*t).*(2-2*exp(-sqrt(-1)*anh.*t));
  
  switch fft_type
    case 'fft'
      %disp('fft')
      S = fftshift(real(fft(S1,n_zp)));
      %S2 = fftshift(real(fft(S2)));
    case 'petersfft'
      %disp('peter''s fft')
      S  = fftshift(real(petersfft(S1,n_interp)));
      %S2 = fftshift(real(petersfft(S2,n_interp)));
    case 'sgrsfft'
      %disp('sgr''s fft')
      S  = fftshift(real(sgrsfft(S1,n_zp)));
      %S2 = fftshift(real(sgrsfft(S2)));
  end
  
 s.c2 = c2(t);
 s.g = g(t);
 s.S1 = S1;
 s.t = t;
 extra = s;
 out = interp1(w,S,w1_in,'pchip');
    return
end

if order==3
    P = zeros(length(w3_in),length(w1_in)); %signle time step
    out = zeros(length(w3_in),length(w1_in),n_t2_array); %array for output
    
  [T1,T3] = meshgrid(t,t);
  for i=1:n_t2_array
    t2 = t2_array(i);
    %letting P1 and P2 be complex seems to work
    %    P1=fft2(exp(sqrt(-1)*w_0.*(-T1+T3)).*exp(-g(T1+phi)+g(t2)-g(T3)-g(T1+phi+t2)-g(t2+T3)+g(T1+phi+t2+T3)).*(2-2.*exp(-sqrt(-1)*anh.*T3)));
    %    P2=fft2(exp(sqrt(-1)*w_0.*(T1+T3)).*exp(-g(T1-phi)-g(t2)-g(T3)+g(T1-phi+t2)+g(t2+T3)-g(T1-phi+t2+T3)).*(2-2.*exp(-sqrt(-1)*anh.*T3)));
    %taking the real part also seems to work
    %    P1=fft2(real(exp(sqrt(-1)*w_0.*(-T1+T3)).*exp(-g(T1+phi)+g(t2)-g(T3)-g(T1+phi+t2)-g(t2+T3)+g(T1+phi+t2+T3)).*(2-2.*exp(-sqrt(-1)*anh.*T3))));
    %    P2=fft2(real(exp(sqrt(-1)*w_0.*(T1+T3)).*exp(-g(T1-phi)-g(t2)-g(T3)+g(T1-phi+t2)+g(t2+T3)-g(T1-phi+t2+T3)).*(2-2.*exp(-sqrt(-1)*anh.*T3))));
    if two_level_system
        %disp('two level system')
        P1=real(exp(-1i*w_0.*(-T1+T3)+1i*phi).*exp(-g(T1)+g(t2)-g(T3)-g(T1+t2)-g(t2+T3)+g(T1+t2+T3)));
        P2=real(exp(-1i*w_0.*(T1+T3)+1i*phi).*exp(-g(T1)-g(t2)-g(T3)+g(T1+t2)+g(t2+T3)-g(T1+t2+T3)));
    else
        %disp('multilevel system')
        %P1=real(exp(1i*w_0.*(-T1+T3)+1i*phi).*exp(-g(T1)+g(t2)-g(T3)-g(T1+t2)-g(t2+T3)+g(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*anh.*T3)));
        %P2=real(exp(1i*w_0.*(T1+T3)+1i*phi).*exp(-g(T1)-g(t2)-g(T3)+g(T1+t2)+g(t2+T3)-g(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*anh.*T3)));
        P1=exp(-g(T1)+g(t2)-g(T3)-g(T1+t2)-g(t2+T3)+g(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*anh.*T3));
        P2=exp(-g(T1)-g(t2)-g(T3)+g(T1+t2)+g(t2+T3)-g(T1+t2+T3)).*(2-2.*exp(-sqrt(-1)*anh.*T3));
    end
    
    if flag_rotating_frame == false
        P1 = exp(1i*w_0.*(-T1+T3)+1i*phi).*P1;
        P2 = exp(1i*w_0.*(T1+T3)+1i*phi).*P2;
    end
    
    if orientational_response
      r = orientationalResponse(tau_R,3,T1,t2,T3);
      P1 = P1.*r;
      P2 = P2.*r;
    end
      
    switch lower(apodization)
      case 'triangular'
        %disp('triangular apodization')
        x = 1:-1/(n_t-1):0;
        window_fxn = zeros(n_t,n_t);
        for j=1:n_t,
          for k = 1:n_t
            window_fxn(j,k) =x(j)*x(k);
          end
        end
        P1 = P1.*window_fxn;
        P2 = P2.*window_fxn;
      case {'gauss','gaussian'}
        %disp('gaussian apodization')
        window_fxn = exp(-(T1.^2+T3.^2)./(2*(n_t*dt/2)^2));
        P1 = P1.*window_fxn;
        P2 = P2.*window_fxn;        
      case 'none'
        %disp('no apodization')
    end
    if noise>0 
      P1 = P1+max(max(real(P1)))*noise.*rand(size(P1));
      P2 = P2+max(max(real(P2)))*noise.*rand(size(P2));
    end
    if flag_plot
     figure(1),clf
     my2dPlot(t,t,real(P1),12)
     figure(2),clf
     my2dPlot(t,t,real(P2),12)
     x = real(P1);
    end
    
    %do fft
    switch fft_type
      case 'fft'
        P1 = fft2(P1,n_zp,n_zp);
        P2 = fft2(P2,n_zp,n_zp);
      case 'petersfft'
        P1=petersfft2(P1,n_interp);
        P2=petersfft2(P2,n_interp);
      case 'sgrsfft'
        P1=sgrsfft2(P1,n_zp);
        P2=sgrsfft2(P2,n_zp);
    end
    
    %change how data packaged
    P=-fftshift(real(fliplr(circshift(P1,[0 -1]))+P2));
    
    %simulate laser bandwidth
    if simulate_bandwidth
      %disp(['2D simulating laser bandwidth ' num2str(bandwidth) ' cm-1 fwhm, for ' num2str(bandwidth_axes-1) ' axis' ]);
      [W1,W3]=meshgrid(w);
      switch bandwidth_axes
        case 3
          BW = exp(-((W1-w_0)+(W3-w_0)).^2 ...
            ./(2*(bandwidth/2.355)^2));
        case 2
          BW = exp(-(W1-w_0).^2 ...
            ./(2*(bandwidth/2.355)^2));
      end
      P = P.*BW;
    else
      %disp('2D not simulating bandwidth');
    end

    [W1,W3] = meshgrid(w,w);
    P = P./max(P(:));
    out(:,:,i) = interp2(W1,W3,P,w1_in,w3_in','*linear');
    
  end %end t2_array loop
end %end if order 3

if order==5
  R = cell(1,n_t2_array);
  for i = 1:n_t2_array
  %set up the time variables
  [T1,T3,T5] = meshgrid(t,t,t);
  t2 = t2_array(i);
  t4 = t4_array(i);

  %cumulant expansion results for response functions
  if two_level_system 
    R1 =exp(1i*w_0.*(T1+T3+T5)+1i*phi).*...
      exp(-g(T1)-g(t2)-g(T3)-g(t4)-g(T5) ...
      +g(T1+t2)+g(t2+T3)+g(T3+t4)+g(t4+T5) ...
      -g(T1+t2+T3)-g(t2+T3+t4)-g(T3+t4+T5) ...
      +g(T1+t2+T3+t4)+g(t2+T3+t4+T5) ...
      -g(T1+t2+T3+t4+T5));
    R2 =exp(1i*w_0.*(-T1+T3+T5)+1i*phi).*...
      exp(-g(T1)+g(t2)-g(T3)-g(t4)-g(T5) ...
      -g(T1+t2)-g(t2+T3)+g(T3+t4)+g(t4+T5) ...
      +g(T1+t2+T3)+g(t2+T3+t4)-g(T3+t4+T5) ...
      -g(T1+t2+T3+t4)-g(t2+T3+t4+T5) ...
      +g(T1+t2+T3+t4+T5));
    R3 =exp(1i*w_0.*(T1-T3+T5)+1i*phi).*...
      exp(-g(T1)+g(t2)-g(T3)+g(t4)-g(T5) ...
      -g(T1+t2)-g(t2+T3)-g(T3+t4)-g(t4+T5) ...
      +g(T1+t2+T3)-g(t2+T3+t4)+g(T3+t4+T5) ...
      +g(T1+t2+T3+t4)+g(t2+T3+t4+T5) ...
      -g(T1+t2+T3+t4+T5));
    R4 =exp(1i*w_0.*(-T1-T3+T5)+1i*phi).*...
      exp(-g(T1)-g(t2)-g(T3)+g(t4)-g(T5) ...
      +g(T1+t2)+g(t2+T3)-g(T3+t4)-g(t4+T5) ...
      -g(T1+t2+T3)+g(t2+T3+t4)+g(T3+t4+T5) ...
      -g(T1+t2+T3+t4)-g(t2+T3+t4+T5) ...
      +g(T1+t2+T3+t4+T5));
  else
    R1 =exp(1i*w_0.*(T1+T3+T5)+1i*phi).*...
      exp(-g(T1)-g(t2)-g(T3)-g(t4)-g(T5) ...
      +g(T1+t2)+g(t2+T3)+g(T3+t4)+g(t4+T5) ...
      -g(T1+t2+T3)-g(t2+T3+t4)-g(T3+t4+T5) ...
      +g(T1+t2+T3+t4)+g(t2+T3+t4+T5) ...
      -g(T1+t2+T3+t4+T5)).* ...
      (4-8*exp(-1i*anh.*(T3+T5)) ...
      -4*exp(-1i*anh.*T5) ...
      +6*exp(-1i*anh.*(T3+2*T5)) ...
      +2*exp(-1i*anh*T3));
    R2 =exp(1i*w_0.*(-T1+T3+T5)+1i*phi).*...
      exp(-g(T1)+g(t2)-g(T3)-g(t4)-g(T5) ...
      -g(T1+t2)-g(t2+T3)+g(T3+t4)+g(t4+T5) ...
      +g(T1+t2+T3)+g(t2+T3+t4)-g(T3+t4+T5) ...
      -g(T1+t2+T3+t4)-g(t2+T3+t4+T5) ...
      +g(T1+t2+T3+t4+T5)).* ...
      (4-8*exp(-1i*anh.*(T3+T5)) ...
      -4*exp(-1i*anh.*T5) ...
      +6*exp(-1i*anh.*(T3+2*T5)) ...
      +2*exp(-1i*anh*T3));
    R3 =exp(1i*w_0.*(T1-T3+T5)+1i*phi).*...
      exp(-g(T1)+g(t2)-g(T3)+g(t4)-g(T5) ...
      -g(T1+t2)-g(t2+T3)-g(T3+t4)-g(t4+T5) ...
      +g(T1+t2+T3)-g(t2+T3+t4)+g(T3+t4+T5) ...
      +g(T1+t2+T3+t4)+g(t2+T3+t4+T5) ...
      -g(T1+t2+T3+t4+T5)).* ...
      (4-8*exp(-1i*anh.*(T5-T3)) ...
      -4*exp(-1i*anh.*T5) ...
      +6*exp(-1i*anh.*(2*T5-T3)) ...
      +2*exp(1i*anh*T3));
    R4 =exp(1i*w_0.*(-T1-T3+T5)+1i*phi).*...
      exp(-g(T1)-g(t2)-g(T3)+g(t4)-g(T5) ...
      +g(T1+t2)+g(t2+T3)-g(T3+t4)-g(t4+T5) ...
      -g(T1+t2+T3)+g(t2+T3+t4)+g(T3+t4+T5) ...
      -g(T1+t2+T3+t4)-g(t2+T3+t4+T5) ...
      +g(T1+t2+T3+t4+T5)).* ...
      (4-8*exp(-1i*anh.*(T5-T3)) ...
      -4*exp(-1i*anh.*T5) ...
      +6*exp(-1i*anh.*(2*T5-T3)) ...
      +2*exp(1i*anh*T3));
  end
  
  % orientational contribution to dephasing
  if orientational_response
    r = orientationalResponse(tau_R,5,T1,t2,T3,t4,T5);
    R1 = R1.*r;
    R2 = R2.*r;
    R3 = R3.*r;
    R4 = R4.*r;
  end
  
  switch lower(apodization)
    case 'triangular'
      %disp('3D triangular apodization')
      x = 1:-1/(n_t-1):0;
      window_fxn = zeros(n_t,n_t,n_t);
      for j=1:n_t,
        for k = 1:n_t
          for l = 1:n_t
            window_fxn(j,k,l) =x(j)*x(k)*x(l);
          end
        end
      end
      R1 = R1.*window_fxn;
      R2 = R2.*window_fxn;
      R3 = R3.*window_fxn;
      R4 = R4.*window_fxn;
    case {'gauss','gaussian'}
      %disp('3D gaussian apodization')
      window_fxn = exp(-(T1.^2+T3.^2)./(2*(n_t*dt/2)^2));
      R1 = R1.*window_fxn;
      R2 = R2.*window_fxn;
      R3 = R3.*window_fxn;
      R4 = R4.*window_fxn;
    case 'none'
      %disp('3D no apodization')
  end
  if noise>0
    R1 = R1+max(max(real(R1)))*noise.*rand(size(R1));
    R2 = R2+max(max(real(R2)))*noise.*rand(size(R2));
    R3 = R1+max(max(real(R3)))*noise.*rand(size(R3));
    R4 = R2+max(max(real(R4)))*noise.*rand(size(R4));
  end

  switch lower(fft_type)
    case 'fft'
      R1 = fftn(R1,[n_zp,n_zp,n_zp]);
      R2 = fftn(R2,[n_zp,n_zp,n_zp]);
      R3 = fftn(R3,[n_zp,n_zp,n_zp]);
      R4 = fftn(R4,[n_zp,n_zp,n_zp]);
    case 'petersfft'
      R1 = petersfft3(R1,n_interp);
      R2 = petersfft3(R2,n_interp);
      R3 = petersfft3(R3,n_interp);
      R4 = petersfft3(R4,n_interp);
    case 'sgrsfft'
      R1 = sgrsfft3(R1,n_zp);
      R2 = sgrsfft3(R2,n_zp);
      R3 = sgrsfft3(R3,n_zp);
      R4 = sgrsfft3(R4,n_zp);      
  end

  %do the flips
  R2 = flipdim(circshift(R2,[0 -1 0]),2);
  R3 = flipdim(circshift(R3,[-1 0 0]),1);
  R4 = flipdim(flipdim(circshift(R4,[-1 -1 0]),1),2);
  R{i} = fftshift(real(R1+R2+R3+R4));
  end %end t2_array loop
  
  if simulate_bandwidth
    %disp(['3D simulating laser bandwidth ' num2str(bandwidth) ...
    %  ' cm-1 fwhm, for ' num2str(bandwidth_axes) ' axes' ]);
    [W1,W3,W5]=meshgrid(w);
    if bandwidth_axes == 3
      BW = exp(-((W1-w_0)+(W3-w_0)+(W5-w_0)).^2./(2*(bandwidth/2.355)^2));
    elseif bandwidth_axes == 2
      BW = exp(-((W1-w_0)+(W3-w_0)).^2./(2*(bandwidth/2.355)^2));
    end      
    R{i} = R{i}.*BW;
  else
    %disp('3D not simulating bandwidth');
  end
  
end %end if order 5



% 
% %
% % Print figures as needed
% %
% if order>=3
%   for i = 1:n_t2_array
%     figure(1000+i)
%     ind = find(w>=range(1)&w<=range(2));
%     my2dPlot(w(ind),w(ind),P{i}(ind,ind),'pumpprobe',false)
%   end
% end
% 
% %%
% if order>=5
% figure(20)
% ind = find(w>=range(1)&w<=range(2));
% for i = 1:n_t2_array,
%   clf
%   my3dPlot2(w(ind),w(ind),w(ind),R{i}(ind,ind,ind),'labels',labels)
%   if n_t2_array>1,pause,end
% end
% end
% %%
% 
% return
% %%
% 
% figure(1)
% subplot(1,2,1)
% plot(t,g(t))
% title('lineshape function g')
% xlabel('t / ps')
% %set(gca,'Xlim',[0 5*tau]);
% subplot(1,2,2)
% plot(t,exp(-g(t)))
% xlabel('t / ps')
% title('exp(-g)')
% %set(gca,'Xlim',[0 5*tau]);
% 
% 
% figure(6),clf
% %offset_3d= R1(1,1,1);
% %R1=R1-offset_3d;
% my3dPlot(w,real(fftshift(R1(2:end,2:end,2:end))))
% 
% figure(7),clf
% %offset_3d = R2(1,1,1);
% %R2=R2-offset_3d;
% my3dPlot(w,real(fftshift(R2(2:end,2:end,2:end))))
% 
% figure(8),clf
% %offset_3d = R3(1,1,1);
% %R3=R3-offset_3d;
% my3dPlot(w,real(fftshift(R3(2:end,2:end,2:end))))
% 
% figure(9),clf
% %offset_3d = R4(1,1,1);
% %R4=R4-offset_3d;
% my3dPlot(w,real(fftshift(R4(2:end,2:end,2:end))))
% 
% 
% %%
% figure(2),clf
% norm_1d = trapz(S)*dw;
% mean_1d = trapz(w.*S)*dw./norm_1d;
% [max_1d,i] = max(S);
% peak_1d = w(i)
% plot(w,S,'-o',[peak_1d peak_1d],[0 max_1d*1.2],'k')
% title('1d spec')
% xlabel('\omega')
% set(gca,'XLim',[-300 300]);
% 
% %%
% figure(10)
% for i = 1:n_t2_array,
%   clf
%   my3dPlot(w,R{i})
%   if n_t2_array>1,pause,end
% end
% 
% %%
% [dummy,ind1] = min((w+100).^2);
% [dummy,ind2] = min((w-100).^2);
% ind = 1:4;
% figure(11),clf
% my3dSlices(w,R,t2_array,t4_array,ind,ind1,ind2)
% 
% figure(12),clf
% my3dProjections(w,R,1)
% 
% %%
% figure(3),clf
% my2dPlot(w,w,fftshift(real(P1(2:end,2:end))),20)
% figure(4),clf
% my2dPlot(w,w,fftshift(real(P2(2:end,2:end))),20)
% figure(5),clf
% ind = find(w>range(1)&w<range(2));
% if isempty(ind), 
%   ind = 1:length(P{1});
% end
% for i = 1:n_t2_array,
%   clf
%   my2dPlot(w(ind),w(ind),P{i}(ind,ind),20);
%   if n_t2_array>1,pause,end
% end
% %%
% pump_probe(ind) = pump_probe(ind)./max(pump_probe(ind));
% for i = 1:n_t2_array,
%   switch projection_type
%     case 'window'
%       projection{i} = sum(P{i}(ind,ind),2);
%     case 'all'
%       projection{i} = sum(P{i},2);
%       projection{i} = projection{i}(ind);
%   end
%   projection{i} = projection{i}./max(projection{i});
% end
% figure(6),
% plot(w(ind),pump_probe(ind),'o',w(ind),projection{1})
% 

