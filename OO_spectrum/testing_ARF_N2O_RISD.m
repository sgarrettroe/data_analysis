%% load dataMatrix
clear all
cd('C:\Users\tyler\Desktop\IR Spectra\Data\numbers\polarization\sepcalave') %set this to data path
load('anBmimV2_03032022sepcalave.mat') %will load in as data, crop_data

ADD_TO_STARTUP

%% HERE FORWARD IS FOR WOBBLING IN A CONE MODELS
% close all
clear dyn

w1 = [2000 2400]; %change these!!!
w3 = [2180 2248];
crop_data(1,1:29) = cropData(Bmimn2oRTv2avesepcal.zzzz,w1,w3);
crop_data(2,1:29) = cropData(Bmimn2oRTv2avesepcal.zzxx,w1,w3);

syms k_u k_d k_hgs10 k_r t


K = [-k_u,       k_r,            k_d,        0;
        0,  -k_r-k_hgs10-k_u,    0,         k_d;
      k_u,       k_hgs10,       -k_d,       k_r;
        0,         k_u,          0,        -k_r-k_d];

[V3, D3] = eig(K);
ev3 = diag(D3);


V3_0 = [1; 0; 0; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

A = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r,t]);
B = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r,t]);
C = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r,t]);
D = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r,t]);

V3_0 = [0; 1; 0; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

E =  matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r,t]);
F = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r,t]);
G = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r,t]);
H = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r,t]);

V3_0 = [0; 0; 1; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

I = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r,t]);
J = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r,t]);
K = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r,t]);
L = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r,t]);

V3_0 = [0; 0; 0; 1];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

M = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r,t]);
N = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r,t]);
O = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r,t]);
P = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r,t]);

%%
AA = @(T1,t2,T3,p) A(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r,t2);
BB = @(T1,t2,T3,p) B(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r,t2);
CC = @(T1,t2,T3,p) C(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r,t2);
DD = @(T1,t2,T3,p) D(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r,t2);
EE = @(T1,t2,T3,p) E(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r,t2);
FF = @(T1,t2,T3,p) F(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r,t2);
GG = @(T1,t2,T3,p) G(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r,t2);
HH = @(T1,t2,T3,p) H(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r,t2);
II = @(T1,t2,T3,p) I(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r,t2);
JJ = @(T1,t2,T3,p) J(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r,t2);
KK = @(T1,t2,T3,p) K(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r,t2);
LL = @(T1,t2,T3,p) L(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r,t2);
MM = @(T1,t2,T3,p) M(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r,t2);
NN = @(T1,t2,T3,p) N(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r,t2);
OO = @(T1,t2,T3,p) O(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r,t2);
PP = @(T1,t2,T3,p) P(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r,t2);


%% ok let's try multiple cones and inertial cone
clear dyn
dyn = aDArrayBnd;%first time call empty, later can call with indices
label = 'kinetics';
fun_array = {AA BB CC DD ... 
             EE FF GG HH ...
             II JJ KK LL...
             MM NN OO PP};
         
ind_array = {[1 4] ,  [],            [13 16],  [],...
            [27 28],  [2 5 3 6],     [25 26],  [14 17 15 18],...
            [19 22],  [],            [7 10],   [],...
            [31 32],  [20 23 21 24], [29 30],  [8 11 9 12]};
dynParams1(1) = fitParamBnd('k_u',    0.007, 0,  1, ''); % from Matt's calc
dynParams1(2) = fitParamBnd('k_hgs10',0.007, 0,  1, ''); %0.008 old value
dynParams1(3) = fitParamBnd('k_r',    0.0113, 0,  1, ''); %vib lifetime (88 ps)
dynParams1(4) = fitParamBnd('k_d',    0.07, 0,  1, ''); %shoulder decay (14 ps)
dyn(1) = aDArrayBnd(label,fun_array,ind_array,dynParams1,'fixed');


%
clear lsfParams lsf

%
scans = [1:length(crop_data)];
clear aRFoptions s;
aRFoptions.dt = 0.2;
aRFoptions.n_t = 96;
aRFoptions.n_zp = 2*aRFoptions.n_t;
aRFoptions.flag_rotating_frame = true;
aRFoptions.w_01_cm =     fitParamBnd('w_01_cm',    2225.7,2220, 2230,'fixed');
aRFoptions.dw_sb_cm =    fitParamBnd('dw_sb_cm',   15,    10,   15,  'fixed');
aRFoptions.mu01sq =      fitParamBnd('mu01sq',     0.082, 1e-5, 1.5, 'fixed');
aRFoptions.mu12sqRatio = fitParamBnd('mu12sqRatio',2   ,  0,    2.5, 'fixed');
aRFoptions.anh_cm =      fitParamBnd('anh_cm',     30,    20,   40,  'fixed');
aRFoptions.phase_deg =   fitParamBnd('phase_deg',  -1.6, -90,   90,  'fixed');
aRFoptions.dyn = dyn;
aRFoptions.dE_cm =       fitParamBnd('dE_cm',      590,   500,  700, 'fixed');
aRFoptions.temperature = fitParamBnd('temperature',295,   290,  505, 'fixed');
%aRFoptions.t2_array = [crop_data(scans).t2]./1000;
%aRFoptions.w1_in = crop_data(1).w1; %2320:0.5:2350;
%aRFoptions.w3_in = (crop_data(1).w3);%2296:2:2350;
aRFoptions.t2_array = [crop_data(1,:).t2]./1000;
aRFoptions.w1_in = crop_data(1,:).w1; 
aRFoptions.w3_in = (crop_data(1,:).w3);
aRFoptions.useParallel=false;  
aRFoptions.nboot = 100;

% wobbling needs lsfOptions and aRFoptions 
% one angles and times are coming from a global fit of the anisotropy data
% at the center of the band
lsfParams(1) =  fitParamBnd('Delta_cm',   6.25,   0.01,  50,     ''); %total linewidth
lsfParams(2) =  fitParamBnd('theta0_deg',   24,   0,     179.999,''); %inertial cone   
lsfParams(3) =  fitParamBnd('tr1',         1.8,   1,      10,    ''); %cone 1 time
lsfParams(4) =  fitParamBnd('theta1_deg',   39,   1,     179.999,''); %cone 1 angle
lsfParams(5) =  fitParamBnd('tr2',          70,   5,     200,    ''); %cone 2
lsfParams(6) =  fitParamBnd('T2',          2.5,   1,     180,    ''); %
lsfParams(7) =  fitParamBnd('ampSSD1',    0.36,   0,       1,    ''); %ssd 1 amp
lsfParams(8) =  fitParamBnd('tauSSD1',    12.7,   1,     100,    ''); %ssd 1 tau
lsfParams(9) =  fitParamBnd('ampSSD2',    0.25,   0,       1,    ''); %ssd 2 amp
lsfParams(10) = fitParamBnd('tauSSD2',     350,   1,     900,    ''); %ssd 2 tau
aRFoptions.pol = 'para';
%aRFoptions.pol = 'perp';
aRFoptions.order = 1;

% the lsf needs params from the aRFoptions also to set up the numerical integration time
lsf = lsfRISDwobblingHL1inertial1cone1diff2SSDNIBnd(lsfParams,'free',aRFoptions);

%Orientation relaxation
label = 'Orientation_rel';
L_l = lsf.L_l;
oPara = @(T1,t2,T3,p) (1/9).*L_l{1}(T1,p).*(1+4/5.*L_l{2}(t2,p)).*L_l{1}(T3,p);
oPerp = @(T1,t2,T3,p) (1/9).*L_l{1}(T1,p).*(1-2/5.*L_l{2}(t2,p)).*L_l{1}(T3,p);
fun_array = {oPara}; 
ind_array = {[1:32]};
dynOther = additionalDynamics(label,fun_array,ind_array);%has no free parameters (gets them from lsf)

%
aRFoptions.damping = lsf;
aRFoptions.dynOther = dynOther;

s(1) = aRFCO2Bnd(aRFoptions);
s(1).dataMatrix = prepareGlobalFitData(crop_data(1,:))/(6e4);

%below should be for perp
aRFoptions.pol = 'perp';
lsf = lsfRISDwobblingHL1inertial1cone1diff2SSDNIBnd(lsfParams,'free',aRFoptions); %have to reinitialize this so that the 'perp' get applied
fun_array = {oPerp}; 
dynOther = additionalDynamics(label,fun_array,ind_array);%has no free parameters (gets them from lsf)
aRFoptions.damping = lsf;
aRFoptions.dynOther = dynOther;
s(2) = aRFCO2Bnd(aRFoptions); %make a second s for perp spectra
s(2).dataMatrix = prepareGlobalFitData(crop_data(2,:))/(1e4);

nscans = zeros(1,length(scans));
for ii = 1:length(scans)
    nscans(ii) = crop_data(scans(ii)).PARAMS.nScans;
end

wt = nScansToWeights(s(1).dataMatrix,nscans);
s(1).weightMatrix = wt;
wt = nScansToWeights(s(2).dataMatrix,nscans);
s(2).weightMatrix = wt;

t=linspace(0,250,100);
p = lsf.copyParamValuesToParamStruct;
figure(29),plot(t,lsf.R.para(t,p)),title('FFCF R para and perp')
hold on
plot(t,lsf.R.perp(t,p))
hold off
figure(30),
for ii = 1:4
    subplot(4,1,ii)
    plot(t,lsf.L_l{ii}(t,p)),set(gca,'Ylim',[0,1]),title(['OCF L_{',num2str(ii),'}'])
end

m = multiMeasurement(s);

%% simulate via multiMeasurement

%parallel
tic
m.s(1) = m.s(1).calcSpectrum(m.s(1).p0);
% matrixPlot2DIR(s(1).simMatrix,s(1).w1_in,s(1).w3_in,s(1).t2_array,[5 6],'height_px',100,'zlimit',1) %sim

%perpendicular
m.s(2) = m.s(2).calcSpectrum(m.s(2).p0);
% matrixPlot2DIR(s(2).simMatrix,s(2).w1_in,s(2).w3_in,s(2).t2_array,[5 6],'height_px',100,'zlimit',1) %sim
toc

%% minimize error?
% m.globalFit %takes waaaaay too long to even start up the calculations?

%% time to get CLS values from the simMatrices?

zzzz_data = struct('w1',num2cell(1:29), 'w3',num2cell(1:29),'R',num2cell(1:29),'t2',num2cell(1:29));

for ii = 1 : 29
zzzz_data(ii).R = m.s(1).simMatrix(:,:,ii);
zzzz_data(ii).w1 = m.s(1).w1_in;
zzzz_data(ii).w3 = m.s(1).w3_in;
zzzz_data(ii).t2 = m.s(1).t2_array(ii).*1000;
zzzz_data(ii).R(isnan(zzzz_data(ii).R))=1e-8; %
end

d_zzzz = cropData(zzzz_data(1:end),w1,w3);
t2_array_zzzz =  [d_zzzz.t2]./1000;

zzxx_data = struct('w1',num2cell(1:29), 'w3',num2cell(1:29),'R',num2cell(1:29),'t2',num2cell(1:29));

for ii = 1 : 29
zzxx_data(ii).R = m.s(2).simMatrix(:,:,ii);
zzxx_data(ii).w1 = m.s(2).w1_in;
zzxx_data(ii).w3 = m.s(2).w3_in;
zzxx_data(ii).t2 = m.s(2).t2_array(ii).*1000;
zzxx_data(ii).R(isnan(zzxx_data(ii).R))=1e-8; %
end
testpara = cropData(zzzz_data,[2180 2250],[2180 2250]);

matrixPlot2DIR(prepareGlobalFitData(testpara),testpara(1).w1,testpara(1).w3,[testpara.t2]./1000,[5 6],'height_px',100,'zlimit',1); %sim

  %% let's try and figure out how the ESA peak intensity changes as a function of time (normalized?)
  
  for ii = 1 : length(t2_array_zzzz)
      I_red_spara(ii) = trapz(zzzz_data(ii).R(6,:))./trapz(zzzz_data(1).R(6,:));
      I_red_sperp(ii) = trapz(zzxx_data(ii).R(6,:))./trapz(zzxx_data(1).R(6,:));
      % 
      I_red_epara(ii) = trapz(crop_data(1,ii).R(6,:))./trapz(crop_data(1,1).R(6,:));
      I_red_eperp(ii) = trapz(crop_data(2,ii).R(6,:))./trapz(crop_data(2,1).R(6,:));
  end
  
  figure(501)
  hold on
  set(gca,'XScale','log')
  plot(t2_array_zzzz,I_red_spara,'ro')
  plot(t2_array_zzzz,I_red_sperp,'rx')
  plot(t2_array_zzzz,I_red_epara,'bo')
  plot(t2_array_zzzz,I_red_eperp,'bx')
      xlabel('t_2 (ps)')
      ylabel('Normalized Intensity')
      legend({'Sim para','Sim perp','Exp para','Exp perp'})
      title('ESA Intensity')
      xlim([0 210])
      ylim([0 1.2])
  hold off
  
  %% let's try and figure out how the GSB/SE peak intensity changes as a function of time (normalized?)
  
  for ii = 1 : length(t2_array_zzzz)
      I_blue_spara_norm(ii) = trapz(zzzz_data(ii).R(20,:))./trapz(zzzz_data(1).R(20,:));
      I_blue_sperp_norm(ii) = trapz(zzxx_data(ii).R(20,:))./trapz(zzxx_data(1).R(20,:));
  
      I_blue_epara_norm(ii) = trapz(crop_data(1,ii).R(20,:))./trapz(crop_data(1,1).R(20,:));
      I_blue_eperp_norm(ii) = trapz(crop_data(2,ii).R(20,:))./trapz(crop_data(2,1).R(20,:));
  end
  
  figure(503)
  hold on
  set(gca,'XScale','log')
  plot(t2_array_zzzz,I_blue_spara_norm,'ro')
  plot(t2_array_zzzz,I_blue_sperp_norm,'rx')
  plot(t2_array_zzzz,I_blue_epara_norm,'bo')
  plot(t2_array_zzzz,I_blue_eperp_norm,'bx')
      xlabel('t_2 (ps)')
      ylabel('Normalized Intensity')
      legend({'Sim para','Sim perp','Exp para','Exp perp'})
      title('GSB/SE Intensity')
      xlim([0 210])
      ylim([0 1.2])
  hold off

  %% check the anisotropy and isotropic decays
p=struct();

p.k_u     = dynParams1(1).value;
p.k_hgs10 = dynParams1(2).value;
p.k_r     = dynParams1(3).value;
p.dE_cm = 590;
p.temperature = 290 ;
p.Delta_cm   = lsfParams(1).value; %total linewidth
p.theta0_deg = lsfParams(2).value; %inertial cone   
p.tr1        = lsfParams(3).value; %cone 1 time
p.theta1_deg = lsfParams(4).value; %cone 1 angle
p.tr2        = lsfParams(5).value; %cone 2
% p.theta2_deg = lsfParams(6).value; %cone 2 angle
% p.tr3        = lsfParams(7).value; %diffusive time or t3

  for ii = 1 : length(t2_array_zzzz)
      I_blue_spara_unnorm(ii) = trapz(zzzz_data(ii).R(20,:));
      I_blue_sperp_unnorm(ii) = trapz(zzxx_data(ii).R(20,:));
      
      I_blue_epara_unnorm(ii) = trapz(crop_data(1,ii).R(20,:));
      I_blue_eperp_unnorm(ii) = trapz(crop_data(2,ii).R(20,:));
  end
   

figure(2001)
hold on
set(gca,'XScale','log')
plot(t2_array_zzzz,(I_blue_epara_unnorm-I_blue_eperp_unnorm)./(I_blue_epara_unnorm+2.*I_blue_eperp_unnorm),'x')
plot(t2_array_zzzz,(I_blue_spara_unnorm-I_blue_sperp_unnorm)./(I_blue_spara_unnorm+2.*I_blue_sperp_unnorm),'o')
plot(bmimtf2nPPn2o.t2,bmimtf2nPPn2o.anisoFit,'^')
plot(t2_array_zzzz,0.4.*L_l{2}(t2_array_zzzz,p))
    xlabel('t_2 (ps)')
    ylabel('anisotropy')
    legend('Exp','Sim','PP','L_2')
hold off


figure(2002)
% semilogx(t2_array_zzzz,(I_blue_epara_norm+2.*I_blue_eperp_norm)./3,'rx')
hold on
set(gca,'XScale','log')
plot(t2_array_zzzz,(I_blue_epara_norm+2.*I_blue_eperp_norm)./3,'x')
plot(bmimtf2nPPn2o.t2,bmimtf2nPPn2o.iso./bmimtf2nPPn2o.iso(1),'^')
plot(t2_array_zzzz,(I_blue_spara_norm+2.*I_blue_sperp_norm)./3,'o')
    xlabel('t_2 (ps)')
    ylabel('isotropic decay')
%     legend('Rpara+2Rperp','L_1')
hold off


%%

HWHM = 6.3;
pfOptions.range1 = [2225.7 - HWHM   2225.7 + HWHM];
pfOptions.range3 = [2180 2248]; 
view2DIRdata(zzzz_data(1),pfOptions.range1,pfOptions.range3);
pfOptions.fitfcn = fittype('- a1.*voigt(w,center,wg,wl) + a2.*voigt(w,center-anh,wg,wl) +c', ...
    'coeff',{'center','a1','a2','wg','wl','anh','c'},'indep',{'w'});
pfOptions.startpoint = [2225.7 1e-1   1e-1    5    5    30   0    ];
pfOptions.lb = [        2210   0      0       0    0    20   0    ];
pfOptions.ub = [        2240   Inf    Inf     100  100  40   1000 ];
pfOptions.zp_factor = 2;
pfOptions.flag_plot = 0;
pfOptions.estimateArea = 1;

d_zzzz = cropData(zzzz_data(1:end),pfOptions.range1,pfOptions.range3);
t2_array_zzzz =  [d_zzzz.t2]./1000;

%%
peakFit_zzzz = extractMaxima(d_zzzz,pfOptions);

%%
maxMatrix_zzzz = CLSMaxMatrix(peakFit_zzzz);
[CLS_zzzz,c2_zzzz,c2_std_zzzz] = fitCLS(d_zzzz,maxMatrix_zzzz,0);

%% ZZXX fitting
HWHM = 6.3;
pfOptions.range1 = [2225.7 - HWHM 2225.7 + HWHM];
pfOptions.range3 = [2180 2248]; 
view2DIRdata(zzzz_data(1),pfOptions.range1,pfOptions.range3);
pfOptions.fitfcn = fittype('- a1.*voigt(w,center,wg,wl) + a2.*voigt(w,center-anh,wg,wl) + c', ...
    'coeff',{'center','a1','a2','wg','wl','anh','c'},'indep',{'w'});
pfOptions.startpoint = [2225.7 1e-1   1e-1    5    5    30   0    ];
pfOptions.lb = [        2220   0      0       0    0    20   0    ];
pfOptions.ub = [        2230   Inf    Inf     100  100  40   1000 ];
pfOptions.zp_factor = 2;
pfOptions.flag_plot = 0;
pfOptions.estimateArea = 1;

d_zzxx = cropData(zzxx_data(1:end),pfOptions.range1,pfOptions.range3);
t2_array_zzxx =  [d_zzxx.t2]./1000;

%%
peakFit_zzxx = extractMaxima(d_zzxx,pfOptions);
%%
maxMatrix_zzxx = CLSMaxMatrix(peakFit_zzxx);
[CLS_zzxx,c2_zzxx,c2_std_zzxx] = fitCLS(d_zzxx,maxMatrix_zzxx,0);

%% Fitting Corr.Funct
% close all

cfOptions(1).startpoint = [0.3  3    0.0 ];
cfOptions(1).lb =         [0    0.1  0    ];
cfOptions(1).ub =         [1    500  1    ];
cfOptions(1).fitfcn = fittype(@(a1,t1,c,x) ...
    a1.*exp(-x./t1) + c, ...
    'coeff', {'a1', 't1', 'c'}, 'indep', {'x'});
cfOptions(1).flag_plot = 0; %0=off 1=on  
% 
cfOptions(2).startpoint = [0.3 0.3 3    100 ];
cfOptions(2).lb =         [0   0   0.1  10  ];
cfOptions(2).ub =         [1   1   50   800 ];
cfOptions(2).fitfcn = fittype(@(a1,a2,t1,t2,x) ...
    a1.*exp(-x./t1) +  a2*exp(-x./t2) , ...
    'coeff', {'a1', 'a2', 't1','t2'}, 'indep', {'x'});
cfOptions(2).flag_plot = 1; %0=off 1=on 

cfOptions(3).startpoint = [0.2 0.2 0.2 3    10  100 ];
cfOptions(3).lb =         [0   0   0   0.1  5   50  ];
cfOptions(3).ub =         [1   1   1   5    50  Inf];
cfOptions(3).fitfcn = fittype(@(a1,a2,a3,t1,t2,t3,x) ...
      a1.*exp(-x./t1) +  a2*exp(-x./t2)+a3.*exp(-x./t3), ...
      'coeff', {'a1', 'a2','a3', 't1','t2','t3'}, 'indep', {'x'});
cfOptions(3).flag_plot = 1; %0=off 1=on

for ii = 1:3
    corr_struct_zzzz = corrFcnFit(t2_array_zzzz,c2_zzzz,c2_std_zzzz,cfOptions(ii));
    CFFit_zzzz(ii) = corr_struct_zzzz;
end
for ii = 1:3
    corr_struct_zzxx = corrFcnFit(t2_array_zzxx,c2_zzxx,c2_std_zzxx,cfOptions(ii));
    CFFit_zzxx(ii) = corr_struct_zzxx;
end

%%

figure(21) %biexponential
  hold on
  set(gca,'XScale','log')
  set(gca,'TickDir','out')
  set(gcf,'color','w')
  box off
  xlabel('\itt\rm_2 (ps)')
  ylabel('CLS')
  % plot(CFFit_zzzz(2).t2(:),(CFFit_zzzz(3).fitresult(CFFit_zzzz(2).t2(:)))); 
  errorbar(CFFit_zzzz(2).t2(:),CFFit_zzzz(2).c2,CFFit_zzzz(2).c2_std,'o','LineWidth',1.5,'Color','k'); 
  % plot(CFFit_zzxx(2).t2(:),(CFFit_zzxx(3).fitresult(CFFit_zzxx(2).t2(:)))); 
  errorbar(CFFit_zzxx(2).t2(:),CFFit_zzxx(2).c2,CFFit_zzxx(2).c2_std,'x','LineWidth',1.5,'Color','k');
  plot(t2_array_zzzz,(lsfParams(7).value.*exp(-t2_array_zzzz./lsfParams(8).value) ...
                    +lsfParams(9).value.*exp(-t2_array_zzzz./lsfParams(10).value)) ...
                    .*(lsf.R.para(t2_array_zzzz,p)),'LineWidth',2)
  plot(t2_array_zzzz,(lsfParams(7).value.*exp(-t2_array_zzzz./lsfParams(8).value) ...
                     +lsfParams(9).value.*exp(-t2_array_zzzz./lsfParams(10).value)) ...
                    .*(lsf.R.perp(t2_array_zzzz,p)),'LineWidth',2)
  % plot(Bmimn2oRTv2avesepcal.analysis.zzzz(2).t2(:),Bmimn2oRTv2avesepcal.analysis.zzzz(2).fitresult(Bmimn2oRTv2avesepcal.analysis.zzzz(2).t2(:)),'b');
  errorbar(Bmimn2oRTv2avesepcal.analysis.zzzz(2).t2(:),Bmimn2oRTv2avesepcal.analysis.zzzz(2).c2(:),Bmimn2oRTv2avesepcal.analysis.zzzz(2).c2_std(:),'ro');
  % plot(Bmimn2oRTv2avesepcal.analysis.zzxx(2).t2(:),Bmimn2oRTv2avesepcal.analysis.zzxx(2).fitresult(Bmimn2oRTv2avesepcal.analysis.zzxx(2).t2(:)),'r');
  errorbar(Bmimn2oRTv2avesepcal.analysis.zzxx(2).t2(:),Bmimn2oRTv2avesepcal.analysis.zzxx(2).c2(:),Bmimn2oRTv2avesepcal.analysis.zzxx(2).c2_std(:),'rx');
%  
