%% Linear functions and parameters

F = 150;               % sampling frequency
dt = 1/F;              % sampling period
T = 1000;              % length of signal
t1 = dt*(0:T-1);       % time vector
N = length(t1);        % number of data points
freq1 = F/N*(0:N-1);
sigma = 0.1;           
homolifetime = 0.05;
omega1 = 2*pi*50;
omega2 = 2*pi*67;
Delta1 = 1;                
Lambda1 = 1/homolifetime;       

% prompt1 = 'Would you like a Gaussian or Lorentzian lineshape?    ';
% shape = input(prompt1,'s');

% if strcmpi(shape,'L')
    g = @(t1) t1/homolifetime;
% elseif strcmpi(shape,'G')
%     g = @(t1) (t1.^2)/(2*(sigma^2));
% else
%     g = @(t1) Delta1^2/Lambda1^2.*(exp(-Lambda1.*t1)-1+Lambda1.*t1);            % as Lambda1 --> inf: t term dominates
%                                                                                 % as Lambda1 --> 0: t^2 term dominates                                                                                                                                                 % starts looking Lorentzian for Lambda < 1
% end

f1 = @(t1) -1i*omega1.*t1;
f2 = @(t1) -1i*omega2.*t1;

lt1 = exp(-f1(t1)).*exp(-1*g(t1));
lt1(1) = lt1(1)/2;
lfft1 = fft(lt1);
lf1 = abs(lfft1);

lt2 = exp(-f2(t1)).*exp(-1*g(t1));
lt2(1) = lt2(1)/2;
lsgrsfft2 = fft(lt2);
lf2 = abs(lsgrsfft2);

%% Third order functions and parameters
Tau2 = 75;
Tau3 = 1000;
t2 = Tau2;
t3 = dt*(0:Tau3-1);
N2 = length(t2);
N3 = length(t3);
freq3 = F/N3*(0:N3-1);
[T1,T3] = meshgrid(t1,t3);

Rt1_1 = exp(-1*f1(T3)+f1(T1)).*exp(-1*g(T1)+g(t2)-g(T3)-g(T1+t2)-g(t2+T3)+g(T1+t2+T3));     % rephasing response function 1
Rf1_1 = real(sgrsfft2(Rt1_1));
Rf1_1(1) = 0.5*Rf1_1(1);
Rf1_1 = fliplr(Rf1_1);                              % mirror the response
Rf1_1 = circshift(Rf1_1,1,2);                       % put zero back on the left

Rt1_2 = exp(-1*f2(T3)+f2(T1)).*exp(-1*g(T1)+g(t2)-g(T3)-g(T1+t2)-g(t2+T3)+g(T1+t2+T3));     % rephasing response function 2
Rf1_2 = real(sgrsfft2(Rt1_2));
Rf1_2(1) = 0.5*Rf1_2(1);
Rf1_2 = fliplr(Rf1_2);                              % mirror the response
Rf1_2 = circshift(Rf1_2,1,2);                       % put zero back on the left

%%

Rt2_1 = exp(-1*f1(T3+T1)).*exp(-1*g(T1)-g(t2)-g(T3)+g(T1+t2)+g(t2+T3)-g(T1+t2+T3));     % non-rephasing response function 1
Rf2_1 = sgrsfft2(Rt2_1);
Rf2_1(1) = 0.5*Rf2_1(1);
Rf2_1 = real(Rf2_1);

Rt2_2 = exp(-1*f2(T3+T1)).*exp(-1*g(T1)-g(t2)-g(T3)+g(T1+t2)+g(t2+T3)-g(T1+t2+T3));     % non-rephasing response function 2
Rf2_2 = sgrsfft2(Rt2_2);
Rf2_2(1) = 0.5*Rf2_2(1);
Rf2_2 = real(Rf2_2);
%% additional parameters for three level system

% prompt2 = 'Calculate for a 2 or a 3 level system?   ';
% levels = input(prompt2,'s');
% 
% if strcmp(levels,'3')
%     omega3 = 2*pi*75;
%     f3 = @(t3) -1i*omega3.*t3;
%     
%     Rt1_2 = -exp(-1*f3(T3)+f1(T1)).*exp(-1*g(T1)+g(t2)-g(T3)-g(T1+t2)-g(t2+T3)+g(T1+t2+T3));     % rephasing response function
%     Rfft1_2 = sgrsfft2(Rt1_2);
%     Rf1_2 = real(Rfft1_2);
%     Rf1_2 = fliplr(Rf1_2);
%     Rf1_2 = circshift(Rf1_2,1,2);
% 
%     Rt2_2 = -exp(-1*f3(T3)-f1(T1)).*exp(-1*g(T1)-g(t2)-g(T3)+g(T1+t2)+g(t2+T3)-g(T1+t2+T3));     % non-rephasing response function
%     Rsgrsfft2_2 = sgrsfft2(Rt2_2);
%     Rf2_2 = real(Rsgrsfft2_2);
% end

%% plot linear response

figure(1)
plot(real(lt1))
xlabel('time')
ylabel('signal')
title('Signal in Time Domain for Two Level System')
  
figure(2)
plot(freq1,real(lf1)+real(lf2))
xlabel('frequency')
ylabel('signal')
title('Signal in Frequency Domain for Two Level System')
    
% figure(3)
% plot(t1,g)
   
%% plot third order response

[F1,F3] = meshgrid(freq1,freq3);
   
% if strcmp(levels,'2')
%      Rsum1 = Rf1_1 + 2*Rf2_1; 
%      Rsum2 = Rf1_2 + 2*Rf2_2; 
% sum non-rephasing and mirror rephasing response to get absorptive correlation
%     rephase = Rf1_1 + Rf1_2;
%     figure(4)
%     contour(F1,F3,rephase)
%     xlabel('frequency 1')
%     ylabel('frequency 3')
%     title('Sum of Rephasing at Two Different Frequencies')
    
%     figure(5)
%     contour(F1,F3,Rf2_1)
%     xlabel('frequency 1')
%     ylabel('frequency 3')
%     title('2D Non-Rephasing Spectrum for Two Level System')
%     
%     figure(6)
%     contour(F1,F3,Rsum1)
%     hold on
%     contour(F1,F3,Rsum2)
%     xlabel('frequency')
%     ylabel('signal')
%     title('Overlap due to Phase Twisting')
%     
%     figure(7)
%     my2dPlot(F1,F3,-(Rsum1+Rsum2))
%     xlabel('frequency')
%     ylabel('signal')
%     title('Phase Twisting due to Imbalance of Pathways')
    %figure(7)
    %contour(T1,T3,real(Rt1_1))
    
    %figure(8)
    %contour(T1,T3,real(Rt2_1))
    
% elseif strcmp(levels,'3')
    Rf1 = Rf1_1; % + Rf1_2;
    Rf2 = Rf2_1; % + Rf2_2;
    Rsum = Rf1 + Rf2;
%     
    figure(4)
    contour(F1,F3,Rf1)
    xlabel('\omega_1/2\pic (cm^{-1})')
    ylabel('\omega_3/2\pic (cm^{-1})')
    title('2D Rephasing Spectrum')
    
    figure(5)
    contour(F1,F3,Rf2)
    xlabel('\omega_1/2\pic (cm^{-1})')
    ylabel('\omega_3/2\pic (cm^{-1})')
    title('2D Non-Rephasing Spectrum')
%     
    figure(6)
    contour(F1,F3,Rsum)
    xlabel('\omega_1/2\pic (cm^{-1})')
    ylabel('\omega_3/2\pic (cm^{-1})')
    title('Absorptive 2D Spectrum')
%     
%     %figure(7)
%     %contour(T1,T3,real(Rt1))
%     
%     %figure(8)
%     %contour(T1,T3,real(Rt2))
% end

%% I(0,0)

R_w0 = Rsum;
contour(freq1,freq3,real(Rsum))
title('I(0,0)')
xlabel('\omega_1/2\pic (cm^{-1})')
ylabel('\omega_3/2\pic (cm^{-1})')
colorbar
%% Get I(pi/2,0)

R_t1 = sgrsifft(Rsum,[],1);
R_t1 = fftshift(R_t1,1);
t1_new = fftTimeAxis(freq1,'time','ps', 'shift', 'on');
w1_new = fftFreqAxis(t1_new, 'time', 'ps', 'shift', 'off');

% zero negative time
R_t1(find(t1_new < 0),:) = 0;

% FFT back to freq
R_w1 = sgrsfft(R_t1,[],1);

contour(w1_new,freq3,real(R_w1))
title('I(\pi/2,0)')
xlabel('\omega_1/2\pic (cm^{-1})')
ylabel('\omega_3/2\pic (cm^{-1})')
colorbar
%% Get I(0,pi/2)

R_t3 = sgrsifft(Rsum,[],2);
R_t3 = fftshift(R_t3,2);
t3_new = fftTimeAxis(freq3,'time','ps', 'shift', 'on');
w3_new = fftFreqAxis(t3_new, 'time', 'ps', 'shift', 'off');

% zero negative time
R_t3(:,find(t3_new < 0)) = 0;

% FFT back to freq
R_w3 = sgrsfft(R_t3,[],2);

contour(w1_new,w3_new,real(R_w3))
title('I(0,\pi/2)')
xlabel('\omega_1/2\pic (cm^{-1})')
ylabel('\omega_3/2\pic (cm^{-1})')
colorbar
%% Get I(pi/2,pi/2)
R_t13 = sgrsifft2(Rsum);
R_t13 = fftshift(R_t13);

% zero negative time
R_t13(find(t1_new < 0),:) = 0;
R_t13(:,find(t3_new < 0)) = 0;

% FFT back to freq
R_w13 = sgrsfft2(R_t13);

contour(w1_new,w3_new,real(R_w13))
title('I(\pi/2,\pi/2)')
xlabel('\omega_1/2\pic (cm^{-1})')
ylabel('\omega_3/2\pic (cm^{-1})')
colorbar
%% Combine to get rephasing spectrum
rephase = R_w0 + R_w13 + 1i*(R_w3 - R_w1);

% plot
contour(w1_new,w3_new,real(rephase))
xlabel('\omega_1/2\pic (cm^{-1})')
ylabel('\omega_3/2\pic (cm^{-1})')
title('Extracted Rephasing Spectrum')
colorbar
%% Combine to get non-rephasing spectrum
nonrephase = R_w0 - R_w13 - 1i*(R_w3 + R_w1);

% plot
contour(w1_new,w3_new,real(nonrephase))
xlabel('\omega_1/2\pic (cm^{-1})')
ylabel('\omega_3/2\pic (cm^{-1})')
title('Extracted Non-Rephasing Spectrum')
colorbar
