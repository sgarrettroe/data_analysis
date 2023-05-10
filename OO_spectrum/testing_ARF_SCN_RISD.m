%% load dataMatrix

ADD_TO_STARTUP
cd('C:\Users\tyler\Desktop\IR Spectra\Scripts\SGR_Shared') %set this to data path
load('LIFE006_Aniso.mat','Aniso_19') %will load in as data, crop_data
load('LIFE006_Aniso.mat','data') 
load('LIFE006_CLS_para.mat','CFFit_CLS_LIFE006_para_19C')
load('LIFE006_CLS_perp.mat','CFFit_CLS_LIFE006_perp_19C')

cd('C:\Users\tyler\Desktop\IR Spectra\Data\clinton')
load('2017_02_15_LIFE006_Tdep.mat')

%% This is for the actual 2dir spectra (Life006)

para19 = Tdep(2).para;
perp19 = Tdep(2).perp;

para19 = sort2DIRdata(para19); %sorts data according to t2
perp19 = sort2DIRdata(perp19);

para19 = para19(7:35); %only taking 200fs to 11ps
perp19 = perp19(7:35);
t2_array = [para19(:).t2]./1000;

w1 = [1990 2130];
w3 = [1990 2130];
crop_data(1,1:29) = cropData(para19,w1,w3);
crop_data(2,1:29) = cropData(perp19,w1,w3);

[W1, W3] = meshgrid(crop_data(1).w1, crop_data(1).w3');
%we have to interpolate the data because they have different w3 axes and
%that messes up the global fitting. 
for ii = 1:29
crop_data(2,ii).R = interp2(crop_data(2,ii).w1, crop_data(2,ii).w3, crop_data(2,ii).R, W1, W3);
crop_data(2,ii).w1 = crop_data(1).w1;
crop_data(2,ii).w3 = crop_data(1).w3;
end 

%% WOBBLING IN ONE CONE from clinton analysis anisotropy values
% close all
clear dyn

% w1 = [1990:3:2100]; %change these!!!
% w3 = [1990:4:2110];
% crop_data(1:28) = cropData(data,w1,w3);
% crop_data(2,1:29) = cropData(Bmimn2oRTv2avesepcal.zzxx,w1,w3);

syms k_u k_d k_hgs10 k_r1 k_r2 t


K = [-k_u,       k_r1,            k_d,        0;
        0,  -k_r1-k_hgs10-k_u,    0,         k_d;
      k_u,       k_hgs10,       -k_d,       k_r2;
        0,         k_u,          0,        -k_r2-k_d];

[V3, D3] = eig(K);
ev3 = diag(D3);


V3_0 = [1; 0; 0; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

A = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
B = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
C = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
D = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

V3_0 = [0; 1; 0; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

E = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
F = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
G = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
H = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

V3_0 = [0; 0; 1; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

I = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
J = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
K = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
L = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

V3_0 = [0; 0; 0; 1];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

M = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
N = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
O = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
P = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

%%
AA = @(T1,t2,T3,p) A(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
BB = @(T1,t2,T3,p) B(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
CC = @(T1,t2,T3,p) C(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
DD = @(T1,t2,T3,p) D(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
EE = @(T1,t2,T3,p) E(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
FF = @(T1,t2,T3,p) F(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
GG = @(T1,t2,T3,p) G(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
HH = @(T1,t2,T3,p) H(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
II = @(T1,t2,T3,p) I(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
JJ = @(T1,t2,T3,p) J(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
KK = @(T1,t2,T3,p) K(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
LL = @(T1,t2,T3,p) L(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
MM = @(T1,t2,T3,p) M(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
NN = @(T1,t2,T3,p) N(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
OO = @(T1,t2,T3,p) O(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
PP = @(T1,t2,T3,p) P(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);


%% ok let's try one inertial, one cone, and diffusion
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
dynParams1(1) = fitParamBnd('k_u',     0.0001, 0,  1, ''); %THESE CANT GO TO ZERO
dynParams1(2) = fitParamBnd('k_hgs10', 0.0001, 0,  1, ''); %THESE CANT GO TO ZERO
dynParams1(3) = fitParamBnd('k_r1',    0.0018, 0,  1, ''); % taken from T1 = 5.7 ps
dynParams1(4) = fitParamBnd('k_r2',    0.0050, 0,  1, ''); % taken from T1 = 2 ps

dyn(1) = aDArrayBnd(label,fun_array,ind_array,dynParams1,'fixed');


%
clear lsfParams lsf

%
% scans = [1:length(crop_data)];
clear aRFoptions s;
aRFoptions.dt = 0.2;
aRFoptions.n_t = 96;
aRFoptions.n_zp = 2*aRFoptions.n_t;
aRFoptions.flag_rotating_frame = true;
aRFoptions.w_01_cm =     fitParamBnd('w_01_cm',    2059,  2020, 2080,'fixed');
aRFoptions.dw_sb_cm =    fitParamBnd('dw_sb_cm',   -11,  -15,  -10,  'fixed');
aRFoptions.mu01sq =      fitParamBnd('mu01sq',     0.08,   1e-5, 1.5, 'fixed');
aRFoptions.mu12sqRatio = fitParamBnd('mu12sqRatio',1.8,   0,    2.5, 'fixed');
aRFoptions.anh_cm =      fitParamBnd('anh_cm',      20,   10,   40,  'fixed');
aRFoptions.phase_deg =   fitParamBnd('phase_deg',    0,  -90,   90,  'fixed');
aRFoptions.dyn = dyn;
aRFoptions.dE_cm =       fitParamBnd('dE_cm',      700,   0,    700, 'fixed');
aRFoptions.temperature = fitParamBnd('temperature',292,   290,  505, 'fixed');
%aRFoptions.t2_array = [crop_data(scans).t2]./1000;
%aRFoptions.w1_in = crop_data(1).w1; %2320:0.5:2350;
%aRFoptions.w3_in = (crop_data(1).w3);%2296:2:2350;
aRFoptions.t2_array = [Aniso_19.t2]./1000;
aRFoptions.w1_in = w1; 
aRFoptions.w3_in = w3; %using Aniso data.w3 requires transpose
aRFoptions.useParallel=false;  
aRFoptions.nboot = 100;

% wobbling needs lsfOptions and aRFoptions 
% one angles and times are coming from a global fit of the anisotropy data
% at the center of the band
lsfParams(1) = fitParamBnd('Delta_cm',      20,   0.01,  50,     ''); %total linewidth
% lsfParams(2) = fitParamBnd('theta0_deg',   15,   0,     180,    ''); %inertial cone   
lsfParams(3) = fitParamBnd('tr1',           22.9,   1,     100,    ''); %cone 1 time 
lsfParams(4) = fitParamBnd('theta1_deg',    56,     1,     179.999,''); %cone 1 angle
% lsfParams(5) = fitParamBnd('tr2',          50,   5,     200,    ''); %cone 2
% lsfParams(6) = fitParamBnd('theta2_deg',   30,   1,     180,    ''); %cone 2 angle
% lsfParams(7) = fitParamBnd('tr3',         100,   50,    500,    ''); %diffusive time or t3
% lsfParams(8) = fitParamBnd('theta3_deg',  90,   1,     180,    ''); %cone 3 angle
% best fit from bootstrapping anisotropy
aRFoptions.pol = 'para';
%aRFoptions.pol = 'perp';
aRFoptions.order = 1;

% the lsf needs params from the aRFoptions also to set up the numerical integration time
lsf = lsfRISDwobbling1coneNIBnd(lsfParams,'free',aRFoptions);

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
% s(1).dataMatrix = prepareGlobalFitData(crop_data(1,:))/(6e4);

%below should be for perp
aRFoptions.pol = 'perp';
lsf = lsfRISDwobbling1coneNIBnd(lsfParams,'free',aRFoptions); %have to reinitialize this so that the 'perp' get applied
fun_array = {oPerp}; 
dynOther = additionalDynamics(label,fun_array,ind_array);%has no free parameters (gets them from lsf)
aRFoptions.damping = lsf;
aRFoptions.dynOther = dynOther;
s(2) = aRFCO2Bnd(aRFoptions); %make a second s for perp spectra
% s(2).dataMatrix = prepareGlobalFitData(crop_data(2,:))/(1.5e4);

% nscans = zeros(1,length(scans));
% for ii = 1:length(scans)
%     nscans(ii) = crop_data(scans(ii)).PARAMS.nScans;
% end

% wt = nScansToWeights(s(1).dataMatrix,nscans);
% s(1).weightMatrix = wt;
% wt = nScansToWeights(s(2).dataMatrix,nscans);
% s(2).weightMatrix = wt;

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
%% chi by eye?
matrixPlot2DIR(m.s(1).simMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit',1); %sim
matrixPlot2DIR(m.s(2).simMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit',1); %sim

%% minimize error?
% m.globalFit %takes waaaaay too long to even start up the calculations?

%% output
% 
% matrixPlot2DIR(m.s(1).simMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit'); %sim
% matrixPlot2DIR(m.s(2).simMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit'); %sim

% matrixPlot2DIR(m.s(1).dataMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit'); %real
% matrixPlot2DIR(m.s(2).dataMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit'); %real
%cool this thing works for once You just need to z-limit the spectra to
%see that they should have the same shape



%% time to get CLS values from the simMatrices?

zzzz_data = struct('w1',num2cell(1:28), 'w3',num2cell(1:28),'R',num2cell(1:28),'t2',num2cell(1:28));

for ii = 1 : 28
zzzz_data(ii).R = m.s(1).simMatrix(:,:,ii);
zzzz_data(ii).w1 = m.s(1).w1_in;
zzzz_data(ii).w3 = m.s(1).w3_in;
zzzz_data(ii).t2 = m.s(1).t2_array(ii).*1000;
end

d_zzzz = cropData(zzzz_data(1:end),w1,w3);
t2_array_zzzz =  [d_zzzz.t2]./1000;

zzxx_data = struct('w1',num2cell(1:28), 'w3',num2cell(1:28),'R',num2cell(1:28),'t2',num2cell(1:28));

for ii = 1 : 28
zzxx_data(ii).R = m.s(2).simMatrix(:,:,ii);
zzxx_data(ii).w1 = m.s(2).w1_in;
zzxx_data(ii).w3 = m.s(2).w3_in;
zzxx_data(ii).t2 = m.s(2).t2_array(ii).*1000;
end
  %% let's try and figure out how the ESA peak intensity changes as a function of time (normalized?)
  
%   for ii = 1 : length(t2_array_zzzz)
%       I_red_spara(ii) = zzzz_data(ii).R(6,22)./zzzz_data(1).R(6,22);
%       I_red_sperp(ii) = zzxx_data(ii).R(6,22)./zzxx_data(1).R(6,22);
%   
%       I_red_epara(ii) = crop_data(1,ii).R(6,21)./crop_data(1,1).R(6,21);
%       I_red_eperp(ii) = crop_data(2,ii).R(6,21)./crop_data(2,1).R(6,21);
%   end
%   
%   figure(501)
%   plot(t2_array_zzzz,I_red_spara,'ro')
%   hold on
%   plot(t2_array_zzzz,I_red_sperp,'rx')
%   plot(t2_array_zzzz,I_red_epara,'bo')
%   plot(t2_array_zzzz,I_red_eperp,'bx')
%       xlabel('t_2 (ps)')
%       ylabel('Normalized Intensity')
%       legend({'Sim para','Sim perp','Exp para','Exp perp'})
%       title('ESA Intensity')
%       xlim([0 210])
%       ylim([0 1.2])
%   hold off
  
  %% let's try and figure out how the GSB/SE peak intensity changes as a function of time (normalized?)
  
  for ii = 1 : length(t2_array_zzzz)
      I_blue_spara_norm(ii) = zzzz_data(ii).R(19,25)./zzzz_data(1).R(19,25);
      I_blue_sperp_norm(ii) = zzxx_data(ii).R(19,25)./zzxx_data(1).R(19,25);
  %center value becomes 70,70 for 1cm-1 resolution
%       I_blue_epara_norm(ii) = crop_data(1,ii).R(20,22)./crop_data(1,1).R(20,22);
%       I_blue_eperp_norm(ii) = crop_data(2,ii).R(20,22)./crop_data(2,1).R(20,22);
  end
  
  figure(101)
  plot(t2_array_zzzz,I_blue_spara_norm,'ro')
  hold on
  plot(t2_array_zzzz,I_blue_sperp_norm,'rx')
%   plot(t2_array_zzzz,I_blue_epara_norm,'bo')
%   plot(t2_array_zzzz,I_blue_eperp_norm,'bx')
      xlabel('t_2 (ps)')
      ylabel('Normalized Intensity')
      legend({'Sim para','Sim perp'})
      title('GSB/SE Intensity')
      xlim([0 12])
      ylim([0 1.2])
  hold off

  %% check the anisotropy and isotropic decays
% p=struct();

% p.k_u     = dynParams1(1).value;
% p.k_hgs10 = dynParams1(2).value;
% p.k_r     = dynParams1(3).value;
% p.dE_cm = 590;
% p.temperature = 290 ;
% p.Delta_cm   = lsfParams(1).value; %total linewidth
% % p.theta0_deg = lsfParams(2).value; %inertial cone   
% p.tr1        = lsfParams(3).value; %cone 1 time
% p.theta1_deg = lsfParams(4).value; %cone 1 angle
% p.tr2        = lsfParams(5).value; %cone 2
% p.theta2_deg = lsfParams(6).value; %cone 2 angle
% p.tr3        = lsfParams(7).value; %diffusive time or t3
% 
%   for ii = 1 : length(t2_array_zzzz)
%       I_blue_spara_unnorm(ii) = zzzz_data(ii).R(20,22);
%       I_blue_sperp_unnorm(ii) = zzxx_data(ii).R(20,22);
%       
%       I_blue_epara_unnorm(ii) = crop_data(1,ii).R(20,22);
%       I_blue_eperp_unnorm(ii) = crop_data(2,ii).R(20,22);
%   end
   
% 
% figure(2001)
% semilogx(t2_array_zzzz,(I_blue_epara_unnorm-I_blue_eperp_unnorm)./(I_blue_epara_unnorm+2.*I_blue_eperp_unnorm),'x')
% hold on
% semilogx(t2_array_zzzz,(I_blue_spara_unnorm-I_blue_sperp_unnorm)./(I_blue_spara_unnorm+2.*I_blue_sperp_unnorm),'o')
% semilogx(t2_array_zzzz,0.4.*L_l{2}(t2_array_zzzz,p))
%     xlabel('t_2 (ps)')
%     ylabel('anisotropy')
%     legend('Exp','Sim','L_2')
% hold off


% figure(2002)
% % semilogx(t2_array_zzzz,(I_blue_epara_norm+2.*I_blue_eperp_norm)./3,'rx')
% 
% semilogx(t2_array_zzzz,(lsf.R.para(t2_array_zzzz,p)+2.*lsf.R.perp(t2_array_zzzz,p))./3,'r')
% hold on
% semilogx(t2_array_zzzz,(I_blue_epara_norm+2.*I_blue_eperp_norm)./3,'ro')
% semilogx(t2_array_zzzz,(I_blue_spara_norm+2.*I_blue_sperp_norm)./3,'bo')
% semilogx(t2_array_zzzz,L_l{1}(t2_array_zzzz,p),'k')
%     xlabel('t_2 (ps)')
%     ylabel('isotropic decay')
% %     legend('Rpara+2Rperp','L_1')
% hold off


%%

HWHM = 9;
pfOptions.range1 = [2040   2082];
pfOptions.range3 = [1990   2100]; 
view2DIRdata(zzzz_data(1),pfOptions.range1,pfOptions.range3);
pfOptions.fitfcn = fittype('- a1.*voigt(w,center,wg,wl) + a2.*voigt(w,center-anh,wg,wl) +c', ...
    'coeff',{'center','a1','a2','wg','wl','anh','c'},'indep',{'w'});
pfOptions.startpoint = [2059   1      1       7    7    15   0    ];
pfOptions.lb = [        2030   0      0       0    0    10   0    ];
pfOptions.ub = [        2090   Inf    Inf     100  100  30   1000 ];
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
HWHM = 9;
pfOptions.range1 = [2040   2082];
pfOptions.range3 = [1990   2100];  
view2DIRdata(zzzz_data(1),pfOptions.range1,pfOptions.range3);
pfOptions.fitfcn = fittype('- a1.*voigt(w,center,wg,wl) + a2.*voigt(w,center-anh,wg,wl) + c', ...
    'coeff',{'center','a1','a2','wg','wl','anh','c'},'indep',{'w'});
pfOptions.startpoint = [2059   1      1       5    5    15   0    ];
pfOptions.lb = [        2030   0      0       0    0    10   0    ];
pfOptions.ub = [        2090   Inf    Inf     100  100  40   1000 ];
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
cfOptions(2).startpoint = [0.3 0.3 3    100 0.0 ];
cfOptions(2).lb =         [0   0   0.1  10  0    ];
cfOptions(2).ub =         [1   1   10   800 1    ];
cfOptions(2).fitfcn = fittype(@(a1,a2,t1,t2,c,x) ...
    a1.*exp(-x./t1) +  a2*exp(-x./t2) + c, ...
    'coeff', {'a1', 'a2', 't1','t2', 'c'}, 'indep', {'x'});
cfOptions(2).flag_plot = 0; %0=off 1=on 

cfOptions(3).startpoint = [0.2 0.2 0.2 3    10  100 ];
cfOptions(3).lb =         [0   0   0   0.1  5   50  ];
cfOptions(3).ub =         [1   1   1   5    50  Inf];
cfOptions(3).fitfcn = fittype(@(a1,a2,a3,t1,t2,t3,x) ...
      a1.*exp(-x./t1) +  a2*exp(-x./t2)+a3.*exp(-x./t3), ...
      'coeff', {'a1', 'a2','a3', 't1','t2','t3'}, 'indep', {'x'});
cfOptions(3).flag_plot = 0; %0=off 1=on

for ii = 1:3
    corr_struct_zzzz = corrFcnFit(t2_array_zzzz,c2_zzzz,c2_std_zzzz,cfOptions(ii));
    CFFit_zzzz(ii) = corr_struct_zzzz;
end
for ii = 1:3
    corr_struct_zzxx = corrFcnFit(t2_array_zzxx,c2_zzxx,c2_std_zzxx,cfOptions(ii));
    CFFit_zzxx(ii) = corr_struct_zzxx;
end

%%

figure(102) %biexponential
  set(gcf,'color','w')
%   plot(CFFit_zzzz(3).t2(:),(CFFit_zzzz(3).fitresult(CFFit_zzzz(3).t2(:)))); 
  hold on
  errorbar(CFFit_zzzz(3).t2(:),CFFit_zzzz(3).c2,CFFit_zzzz(3).c2_std,'bo'); 
  
%   plot(CFFit_zzxx(3).t2(:),(CFFit_zzxx(3).fitresult(CFFit_zzxx(3).t2(:)))); 
  errorbar(CFFit_zzxx(3).t2(:),CFFit_zzxx(3).c2,CFFit_zzxx(3).c2_std,'bx');
%   plot(t2_array_zzzz,lsf.R.para(t2_array_zzzz,p))
%   plot(t2_array_zzzz,lsf.R.perp(t2_array_zzzz,p))
plot(CFFit_CLS_LIFE006_para_19C.t2(:),CFFit_CLS_LIFE006_para_19C.fitresult(CFFit_CLS_LIFE006_para_19C.t2(:)),'k-');
%   errorbar(CFFit_CLS_LIFE006_para_19C.t2(:),CFFit_CLS_LIFE006_para_19C.c2(:),CFFit_CLS_LIFE006_para_19C.c2_std(:),'k');
%   
  plot(CFFit_CLS_LIFE006_perp_19C.t2(:),CFFit_CLS_LIFE006_perp_19C.fitresult(CFFit_CLS_LIFE006_perp_19C.t2(:)),'k--');
%   errorbar(CFFit_CLS_LIFE006_perp_19C.t2(:),CFFit_CLS_LIFE006_perp_19C.c2(:),CFFit_CLS_LIFE006_perp_19C.c2_std(:),'k');
%  
uh = (0.83.*exp(-(CFFit_zzzz(3).t2)./16.5)).*CFFit_zzzz(3).c2;
uhhh = (0.83.*exp(-(CFFit_zzzz(3).t2)./16.5)).*CFFit_zzxx(3).c2;
plot(CFFit_zzzz(3).t2,uh,'ro')
plot(CFFit_zzxx(3).t2,uhhh,'rx')

xlabel('t_2 (ps)')
ylabel('c_2')
hold off
%% do some anisotropy

% for ii = 1:length(t2_array_zzzz)
% paramax(ii) = Bmimn2oRTv2avesepcal.g.para(ii).R(6,9);
% perpmax(ii) = Bmimn2oRTv2avesepcal.g.perp(ii).R(6,9); %all for real data
% 
% paracenter(ii) = Bmimn2oRTv2avesepcal.g.para(ii).R(5,7);
% perpcenter(ii) = Bmimn2oRTv2avesepcal.g.perp(ii).R(5,7);
% end
% 
% anisomax = (paramax - perpmax)./(paramax + 2.*perpmax);
% anisopeak = (paracenter - perpcenter)./(paracenter + 2.*perpcenter);
% t2 = t2_array_zzzz;

%% anisotropy of sim data
w1 = [2058 2061];
w3 = [2058 2061];
simpara = (prepareGlobalFitData(cropData(zzzz_data(1:end),w1,w3)));
simperp = (prepareGlobalFitData(cropData(zzxx_data(1:end),w1,w3)));

simaniso = (simpara - simperp)./(simpara + 2.*simperp);
simanisopeak = squeeze(simaniso(1,1,:))';
% figure(300)
%         [C,h]=contourf(cropData(zzzz_data(1),w1,w3).w1,cropData(zzzz_data(1),w1,w3).w3,simaniso(:,:,1),(0:0.02:0.7));
%         set(gcf,'PaperUnits','centimeters','PaperPosition',[5.0 5.0 7.5 5.0]);
%         title(sprintf('%0.1f ps',t2(1)));
%         xlabel('\omega_1/2\pic (cm^{-1})');
%         ylabel('\omega_3/2\pic (cm^{-1})');
%         clabel(C,h);

%%

figure(103)
set(gcf,'color','w')
hold on
plot([Aniso_19.t2]./1000,[Aniso_19.area],'x')
plot(t2_array_zzzz,simanisopeak,'o')
plot(t2_array_zzzz,0.4.*lsf.L_l{2}(t2_array_zzzz,p),'k')

xlim([0 12])
ylim([0 0.5])
xlabel('t_2 (ps)')
ylabel('anisotropy')
hold off
% 
% lsfParams(1) = fitParamBnd('Delta_cm',     6,   0.01,  50, ''); %total linewidth
% lsfParams(2) = fitParamBnd('theta0_deg',   6,   0,     180,''); %inertial cone   
% lsfParams(3) = fitParamBnd('tr1',          2,   1,     100,''); %cone 1 time
% lsfParams(4) = fitParamBnd('theta1_deg',  27,   1,     180,''); %cone 1 angle
% lsfParams(5) = fitParamBnd('tr2',         12,   5,     200,''); %cone 2
% lsfParams(6) = fitParamBnd('theta2_deg',  45,   1,     180,''); %cone 2 angle
% lsfParams(7) = fitParamBnd('tr3',        120,   50,    500,''); %diffusive time or t3
% lsfParams(8) = fitParamBnd('theta3_deg',90,   1,     180,''); %cone 3 angle
% lsf = lsfRISDwobbling1inertial1cone1diffNIBnd(lsfParams,'free',aRFoptions);
% p = lsf.copyParamValuesToParamStruct;
% plot(t2,0.4.*lsf.L_l{2}(t2,p))
% 
% 
% figure(151)
% semilogx(t2,anisopeak,'x')
% hold on
% semilogx(t2,simanisopeak,'o')
% semilogx(t2,0.4.*lsf.L_l{2}(t2,p))

%% WHAT IF FOR INERTIAL, ONE CONE, AND DIFFUSION
 




%% HERE FORWARD IS FOR WOBBLING IN A SINGLE CONE USING CAJ NUMBERS
% close all
clear dyn

w1 = [1990:3:2100]; %change these!!!
w3 = [1990:4:2110];

% crop_data(1,1:5) = cropData(Bmimn2oRTv2avesepcal.zzzz([1 13 21 23 27]),w1,w3);
% crop_data(2,1:5) = cropData(Bmimn2oRTv2avesepcal.zzxx([1 13 21 23 27]),w1,w3);

syms k_u k_d k_hgs10 k_r1 k_r2 t


K = [-k_u,       k_r1,            k_d,        0;
        0,  -k_r1-k_hgs10-k_u,    0,         k_d;
      k_u,       k_hgs10,       -k_d,       k_r2;
        0,         k_u,          0,        -k_r2-k_d];

[V3, D3] = eig(K);
ev3 = diag(D3);


V3_0 = [1; 0; 0; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

A = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
B = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
C = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
D = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

V3_0 = [0; 1; 0; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

E = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
F = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
G = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
H = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

V3_0 = [0; 0; 1; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

I = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
J = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
K = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
L = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

V3_0 = [0; 0; 0; 1];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

M = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
N = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
O = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
P = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

%%
AA = @(T1,t2,T3,p) A(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
BB = @(T1,t2,T3,p) B(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
CC = @(T1,t2,T3,p) C(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
DD = @(T1,t2,T3,p) D(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
EE = @(T1,t2,T3,p) E(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
FF = @(T1,t2,T3,p) F(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
GG = @(T1,t2,T3,p) G(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
HH = @(T1,t2,T3,p) H(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
II = @(T1,t2,T3,p) I(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
JJ = @(T1,t2,T3,p) J(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
KK = @(T1,t2,T3,p) K(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
LL = @(T1,t2,T3,p) L(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
MM = @(T1,t2,T3,p) M(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
NN = @(T1,t2,T3,p) N(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
OO = @(T1,t2,T3,p) O(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
PP = @(T1,t2,T3,p) P(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);


%% ok let's try one inertial, one cone, and diffusion
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
dynParams1(1) = fitParamBnd('k_u',     0.0001, 0,  1, ''); %THESE CANT GO TO ZERO
dynParams1(2) = fitParamBnd('k_hgs10', 0.0001, 0,  1, ''); %THESE CANT GO TO ZERO
dynParams1(3) = fitParamBnd('k_r1',    0.0018, 0,  1, ''); % taken from T1 = 5.7 ps
dynParams1(4) = fitParamBnd('k_r2',    0.0050, 0,  1, ''); % taken from T1 = 2 ps

dyn(1) = aDArrayBnd(label,fun_array,ind_array,dynParams1,'fixed');

         
clear aRFoptions s;
aRFoptions.dt = 0.2;
aRFoptions.n_t = 96;
aRFoptions.n_zp = 2*aRFoptions.n_t;
aRFoptions.flag_rotating_frame = true;
aRFoptions.w_01_cm =     fitParamBnd('w_01_cm',   2059,   2020, 2080,'fixed');
aRFoptions.dw_sb_cm =    fitParamBnd('dw_sb_cm',   -11,  -15,  -10,  'fixed');
aRFoptions.mu01sq =      fitParamBnd('mu01sq',    0.08,   1e-5, 1.5, 'fixed');
aRFoptions.mu12sqRatio = fitParamBnd('mu12sqRatio',1.8,   0,    2.5, 'fixed');
aRFoptions.anh_cm =      fitParamBnd('anh_cm',      20,   10,   40,  'fixed');
aRFoptions.phase_deg =   fitParamBnd('phase_deg',    0,  -90,   90,  'fixed');
aRFoptions.dyn = dyn;
aRFoptions.dE_cm =       fitParamBnd('dE_cm',      700,   0,   1000, 'fixed'); %amp2 in aRFEAN
aRFoptions.temperature = fitParamBnd('temperature',292,   290,  505, 'fixed');
%aRFoptions.t2_array = [crop_data(scans).t2]./1000;
%aRFoptions.w1_in = crop_data(1).w1; %2320:0.5:2350;
%aRFoptions.w3_in = (crop_data(1).w3);%2296:2:2350;
aRFoptions.t2_array = [Aniso_19.t2]./1000;
aRFoptions.w1_in = w1; 
aRFoptions.w3_in = w3; %using Aniso data.w3 requires transpose
aRFoptions.useParallel=false;  
aRFoptions.nboot = 100;

% wobbling needs lsfOptions and aRFoptions 
% one angles and times are coming from a global fit of the anisotropy data
% at the center of the band
lsfParams(1) = fitParamBnd('Delta_cm',      20,     0.01,  50,     ''); %total linewidth
% lsfParams(2) = fitParamBnd('theta0_deg',     5,     0,     180,    ''); %inertial cone   
lsfParams(3) = fitParamBnd('tr1',          2.6,     1,     100,    ''); %cone 1 time 
lsfParams(4) = fitParamBnd('theta1_deg',   108,     1,     179.999,''); %cone 1 angle
% lsfParams(5) = fitParamBnd('tr2',           70,     5,     200,    ''); %cone 2
% lsfParams(6) = fitParamBnd('theta2_deg',   30,   1,     180,    ''); %cone 2 angle
% lsfParams(7) = fitParamBnd('tr3',         100,   50,    500,    ''); %diffusive time or t3
% lsfParams(8) = fitParamBnd('theta3_deg',  90,   1,     180,    ''); %cone 3 angle
% best fit from clinton paper 1
aRFoptions.pol = 'para';
%aRFoptions.pol = 'perp';
aRFoptions.order = 1;

% the lsf needs params from the aRFoptions also to set up the numerical integration time
lsf = lsfRISDwobbling1coneNIBnd(lsfParams,'free',aRFoptions);

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
% s(1).dataMatrix = prepareGlobalFitData(crop_data(1,:))/(6e4);

%below should be for perp
aRFoptions.pol = 'perp';
lsf = lsfRISDwobbling1coneNIBnd(lsfParams,'free',aRFoptions); %have to reinitialize this so that the 'perp' get applied
fun_array = {oPerp}; 
dynOther = additionalDynamics(label,fun_array,ind_array);%has no free parameters (gets them from lsf)
aRFoptions.damping = lsf;
aRFoptions.dynOther = dynOther;
s(2) = aRFCO2Bnd(aRFoptions); %make a second s for perp spectra
% s(2).dataMatrix = prepareGlobalFitData(crop_data(2,:))/(1.5e4);

t=linspace(0,14,100);
p = lsf.copyParamValuesToParamStruct;
figure(31),plot(t,lsf.R.para(t,p)),title('FFCF R para and perp')
hold on
plot(t,lsf.R.perp(t,p))
hold off
figure(32),
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
%% chi by eye?
matrixPlot2DIR(m.s(1).simMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit',1); %sim
matrixPlot2DIR(m.s(2).simMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit',1) %sim

%% minimize error?
% m.globalFit %takes waaaaay too long to even start up the calculations?

%% output
% % 
% matrixPlot2DIR(m.s(1).simMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit'); %sim
% matrixPlot2DIR(m.s(2).simMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit'); %sim

% matrixPlot2DIR(m.s(1).dataMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit'); %real
% matrixPlot2DIR(m.s(2).dataMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit'); %real
%cool this thing works for once You just need to z-limit the spectra to
%see that they should have the same shape



%% time to get CLS values from the simMatrices?

zzzz_data = struct('w1',num2cell(1:28), 'w3',num2cell(1:28),'R',num2cell(1:28),'t2',num2cell(1:28));

for ii = 1 : 28
zzzz_data(ii).R = m.s(1).simMatrix(:,:,ii);
zzzz_data(ii).w1 = m.s(1).w1_in;
zzzz_data(ii).w3 = m.s(1).w3_in;
zzzz_data(ii).t2 = m.s(1).t2_array(ii).*1000;
end

d_zzzz = cropData(zzzz_data(1:end),w1,w3);
t2_array_zzzz =  [d_zzzz.t2]./1000;

zzxx_data = struct('w1',num2cell(1:28), 'w3',num2cell(1:28),'R',num2cell(1:28),'t2',num2cell(1:28));

for ii = 1 : 28
zzxx_data(ii).R = m.s(2).simMatrix(:,:,ii);
zzxx_data(ii).w1 = m.s(2).w1_in;
zzxx_data(ii).w3 = m.s(2).w3_in;
zzxx_data(ii).t2 = m.s(2).t2_array(ii).*1000;
end
  %% let's try and figure out how the ESA peak intensity changes as a function of time (normalized?)
  
%   for ii = 1 : length(t2_array_zzzz)
%       I_red_spara(ii) = zzzz_data(ii).R(6,22)./zzzz_data(1).R(6,22);
%       I_red_sperp(ii) = zzxx_data(ii).R(6,22)./zzxx_data(1).R(6,22);
%   
%       I_red_epara(ii) = crop_data(1,ii).R(6,21)./crop_data(1,1).R(6,21);
%       I_red_eperp(ii) = crop_data(2,ii).R(6,21)./crop_data(2,1).R(6,21);
%   end
%   
%   figure(501)
%   plot(t2_array_zzzz,I_red_spara,'ro')
%   hold on
%   plot(t2_array_zzzz,I_red_sperp,'rx')
%   plot(t2_array_zzzz,I_red_epara,'bo')
%   plot(t2_array_zzzz,I_red_eperp,'bx')
%       xlabel('t_2 (ps)')
%       ylabel('Normalized Intensity')
%       legend({'Sim para','Sim perp','Exp para','Exp perp'})
%       title('ESA Intensity')
%       xlim([0 210])
%       ylim([0 1.2])
%   hold off
  
  %% let's try and figure out how the GSB/SE peak intensity changes as a function of time (normalized?)
  
  for ii = 1 : length(t2_array_zzzz)
      I_blue_spara_norm(ii) = zzzz_data(ii).R(20,26)./zzzz_data(1).R(20,26);
      I_blue_sperp_norm(ii) = zzxx_data(ii).R(20,26)./zzxx_data(1).R(20,26);
  %center value becomes 70,70 for 1cm-1 resolution
%       I_blue_epara_norm(ii) = crop_data(1,ii).R(20,22)./crop_data(1,1).R(20,22);
%       I_blue_eperp_norm(ii) = crop_data(2,ii).R(20,22)./crop_data(2,1).R(20,22);
  end
  
  figure(201)
  plot(t2_array_zzzz,I_blue_spara_norm,'ro')
  hold on
  plot(t2_array_zzzz,I_blue_sperp_norm,'rx')
%   plot(t2_array_zzzz,I_blue_epara_norm,'bo')
%   plot(t2_array_zzzz,I_blue_eperp_norm,'bx')
      xlabel('t_2 (ps)')
      ylabel('Normalized Intensity')
      legend({'Sim para','Sim perp'})
      title('GSB/SE Intensity')
      xlim([0 12])
      ylim([0 1.2])
  hold off

  %% check the anisotropy and isotropic decays
% p=struct();

% p.k_u     = dynParams1(1).value;
% p.k_hgs10 = dynParams1(2).value;
% p.k_r     = dynParams1(3).value;
% p.dE_cm = 590;
% p.temperature = 290 ;
% p.Delta_cm   = lsfParams(1).value; %total linewidth
% % p.theta0_deg = lsfParams(2).value; %inertial cone   
% p.tr1        = lsfParams(3).value; %cone 1 time
% p.theta1_deg = lsfParams(4).value; %cone 1 angle
% p.tr2        = lsfParams(5).value; %cone 2
% p.theta2_deg = lsfParams(6).value; %cone 2 angle
% p.tr3        = lsfParams(7).value; %diffusive time or t3
% 
%   for ii = 1 : length(t2_array_zzzz)
%       I_blue_spara_unnorm(ii) = zzzz_data(ii).R(20,22);
%       I_blue_sperp_unnorm(ii) = zzxx_data(ii).R(20,22);
%       
%       I_blue_epara_unnorm(ii) = crop_data(1,ii).R(20,22);
%       I_blue_eperp_unnorm(ii) = crop_data(2,ii).R(20,22);
%   end
   
% 
% figure(2001)
% semilogx(t2_array_zzzz,(I_blue_epara_unnorm-I_blue_eperp_unnorm)./(I_blue_epara_unnorm+2.*I_blue_eperp_unnorm),'x')
% hold on
% semilogx(t2_array_zzzz,(I_blue_spara_unnorm-I_blue_sperp_unnorm)./(I_blue_spara_unnorm+2.*I_blue_sperp_unnorm),'o')
% semilogx(t2_array_zzzz,0.4.*L_l{2}(t2_array_zzzz,p))
%     xlabel('t_2 (ps)')
%     ylabel('anisotropy')
%     legend('Exp','Sim','L_2')
% hold off


% figure(2002)
% % semilogx(t2_array_zzzz,(I_blue_epara_norm+2.*I_blue_eperp_norm)./3,'rx')
% 
% semilogx(t2_array_zzzz,(lsf.R.para(t2_array_zzzz,p)+2.*lsf.R.perp(t2_array_zzzz,p))./3,'r')
% hold on
% semilogx(t2_array_zzzz,(I_blue_epara_norm+2.*I_blue_eperp_norm)./3,'ro')
% semilogx(t2_array_zzzz,(I_blue_spara_norm+2.*I_blue_sperp_norm)./3,'bo')
% semilogx(t2_array_zzzz,L_l{1}(t2_array_zzzz,p),'k')
%     xlabel('t_2 (ps)')
%     ylabel('isotropic decay')
% %     legend('Rpara+2Rperp','L_1')
% hold off


%%

% HWHM = 11;
pfOptions.range1 = [2040   2082]; %paper range 
pfOptions.range3 = [1990   2100]; 
view2DIRdata(zzzz_data(1),pfOptions.range1,pfOptions.range3);
pfOptions.fitfcn = fittype('- a1.*voigt(w,center,wg,wl) + a2.*voigt(w,center-anh,wg,wl) +c', ...
    'coeff',{'center','a1','a2','wg','wl','anh','c'},'indep',{'w'});
pfOptions.startpoint = [2059   0.1    0.1     7    7    15   0    ];
pfOptions.lb = [        2030   0      0       0    0    5    0    ];
pfOptions.ub = [        2080   Inf    Inf     100  100  50   1000 ];
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
HWHM = 11;
pfOptions.range1 = [2040   2082];
pfOptions.range3 = [1990   2110];  
view2DIRdata(zzzz_data(1),pfOptions.range1,pfOptions.range3);
pfOptions.fitfcn = fittype('- a1.*voigt(w,center,wg,wl) + a2.*voigt(w,center-anh,wg,wl) + c', ...
    'coeff',{'center','a1','a2','wg','wl','anh','c'},'indep',{'w'});
pfOptions.startpoint = [2059   1      1       5    5    15   0    ];
pfOptions.lb = [        2020   0      0       0    0    10   0    ];
pfOptions.ub = [        2090   Inf    Inf     100  100  40   1000 ];
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
cfOptions(2).startpoint = [0.3 0.3 3    100 0.0 ];
cfOptions(2).lb =         [0   0   0.1  10  0    ];
cfOptions(2).ub =         [1   1   10   800 1    ];
cfOptions(2).fitfcn = fittype(@(a1,a2,t1,t2,c,x) ...
    a1.*exp(-x./t1) +  a2*exp(-x./t2) + c, ...
    'coeff', {'a1', 'a2', 't1','t2', 'c'}, 'indep', {'x'});
cfOptions(2).flag_plot = 0; %0=off 1=on 

cfOptions(3).startpoint = [0.2 0.2 0.2 3    10  100 ];
cfOptions(3).lb =         [0   0   0   0.1  5   50  ];
cfOptions(3).ub =         [1   1   1   5    50  Inf];
cfOptions(3).fitfcn = fittype(@(a1,a2,a3,t1,t2,t3,x) ...
      a1.*exp(-x./t1) +  a2*exp(-x./t2)+a3.*exp(-x./t3), ...
      'coeff', {'a1', 'a2','a3', 't1','t2','t3'}, 'indep', {'x'});
cfOptions(3).flag_plot = 0; %0=off 1=on

for ii = 1:3
    corr_struct_zzzz = corrFcnFit(t2_array_zzzz,c2_zzzz,c2_std_zzzz,cfOptions(ii));
    CFFit_zzzz(ii) = corr_struct_zzzz;
end
for ii = 1:3
    corr_struct_zzxx = corrFcnFit(t2_array_zzxx,c2_zzxx,c2_std_zzxx,cfOptions(ii));
    CFFit_zzxx(ii) = corr_struct_zzxx;
end

%%

figure(202) %biexponential
  set(gcf,'color','w')
%   plot(CFFit_zzzz(3).t2(:),(CFFit_zzzz(3).fitresult(CFFit_zzzz(3).t2(:)))); 
  hold on
  errorbar(CFFit_zzzz(3).t2(:),CFFit_zzzz(3).c2,CFFit_zzzz(3).c2_std,'bo'); 
  
%   plot(CFFit_zzxx(3).t2(:),(CFFit_zzxx(3).fitresult(CFFit_zzxx(3).t2(:)))); 
  errorbar(CFFit_zzxx(3).t2(:),CFFit_zzxx(3).c2,CFFit_zzxx(3).c2_std,'bo');
%   plot(t2_array_zzzz,lsf.R.para(t2_array_zzzz,p))
%   plot(t2_array_zzzz,lsf.R.perp(t2_array_zzzz,p))
plot(CFFit_CLS_LIFE006_para_19C.t2(:),CFFit_CLS_LIFE006_para_19C.fitresult(CFFit_CLS_LIFE006_para_19C.t2(:)),'k-');
%   errorbar(Bmimn2oRTv2avesepcal.analysis.zzzz(3).t2(:),Bmimn2oRTv2avesepcal.analysis.zzzz(3).c2(:),Bmimn2oRTv2avesepcal.analysis.zzzz(3).c2_std(:),'bo');
%   
  plot(CFFit_CLS_LIFE006_perp_19C.t2(:),CFFit_CLS_LIFE006_perp_19C.fitresult(CFFit_CLS_LIFE006_perp_19C.t2(:)),'k--');
%   errorbar(Bmimn2oRTv2avesepcal.analysis.zzxx(3).t2(:),Bmimn2oRTv2avesepcal.analysis.zzxx(3).c2(:),Bmimn2oRTv2avesepcal.analysis.zzxx(3).c2_std(:),'ro');
%  
% uh = 0.83.*CFFit_zzzz(3).c2;
% uhhh = 0.83.*CFFit_zzxx(3).c2;

% 
uh = (0.86.*exp(-(CFFit_zzzz(3).t2)./25)).*CFFit_zzzz(3).c2;
uhhh = (0.86.*exp(-(CFFit_zzzz(3).t2)./25)).*CFFit_zzxx(3).c2; %using SSD term from paper

% uh = (0.2.*exp(-(CFFit_zzzz(3).t2)./5)+0.63.*exp(-(CFFit_zzzz(3).t2)./300)).*CFFit_zzzz(3).c2;
% uhhh = (0.2.*exp(-(CFFit_zzzz(3).t2)./5)+0.63.*exp(-(CFFit_zzzz(3).t2)./300)).*CFFit_zzxx(3).c2;
plot(CFFit_zzzz(3).t2,uh,'ro')
plot(CFFit_zzxx(3).t2,uhhh,'rx')
hold off
%% do some anisotropy

% for ii = 1:length(t2_array_zzzz)
% paramax(ii) = Bmimn2oRTv2avesepcal.g.para(ii).R(6,9);
% perpmax(ii) = Bmimn2oRTv2avesepcal.g.perp(ii).R(6,9); %all for real data
% 
% paracenter(ii) = Bmimn2oRTv2avesepcal.g.para(ii).R(5,7);
% perpcenter(ii) = Bmimn2oRTv2avesepcal.g.perp(ii).R(5,7);
% end
% 
% anisomax = (paramax - perpmax)./(paramax + 2.*perpmax);
% anisopeak = (paracenter - perpcenter)./(paracenter + 2.*perpcenter);
% t2 = t2_array_zzzz;

%% anisotropy of sim data
w1 = [2058 2060];
w3 = [2058 2060];
simpara = (prepareGlobalFitData(cropData(zzzz_data(1:end),w1,w3)));
simperp = (prepareGlobalFitData(cropData(zzxx_data(1:end),w1,w3)));

simaniso = (simpara - simperp)./(simpara + 2.*simperp);
simanisopeak = squeeze(simaniso(1,1,:))';
% figure(300)
%         [C,h]=contourf(cropData(zzzz_data(1),w1,w3).w1,cropData(zzzz_data(1),w1,w3).w3,simaniso(:,:,1),(0:0.02:0.7));
%         set(gcf,'PaperUnits','centimeters','PaperPosition',[5.0 5.0 7.5 5.0]);
%         title(sprintf('%0.1f ps',t2(1)));
%         xlabel('\omega_1/2\pic (cm^{-1})');
%         ylabel('\omega_3/2\pic (cm^{-1})');
%         clabel(C,h);

%%

figure(203)
plot([Aniso_19.t2]./1000,[Aniso_19.area],'x')
hold on
plot(t2_array_zzzz,simanisopeak,'o')
plot(t2_array_zzzz,0.4.*lsf.L_l{2}(t2_array_zzzz,p),'k-')
xlim([0 12])
ylim([0 0.5])
% figure(300)
%         [C,h]=contourf(cropData(zzzz_data(1),w1,w3).w1,cropData(zzzz_data(1),w1,w3).w3,simaniso(:,:,1),(0:0.02:0.7));
%         set(gcf,'PaperUnits','centimeters','PaperPosition',[5.0 5.0 7.5 5.0]);
%         title(sprintf('%0.1f ps',t2(1)));
%         xlabel('\omega_1/2\pic (cm^{-1})');
%         ylabel('\omega_3/2\pic (cm^{-1})');
%         clabel(C,h);

%%






%% HERE FORWARD IS FOR Wobbling in one inertial, and one cone + SSD
% close all
clear dyn

w1 = [1990:3:2100]; %change these!!!
w3 = [1990:4:2110];

% crop_data(1,1:5) = cropData(Bmimn2oRTv2avesepcal.zzzz([1 13 21 23 27]),w1,w3);
% crop_data(2,1:5) = cropData(Bmimn2oRTv2avesepcal.zzxx([1 13 21 23 27]),w1,w3);

syms k_u k_d k_hgs10 k_r1 k_r2 t


K = [-k_u,       k_r1,            k_d,        0;
        0,  -k_r1-k_hgs10-k_u,    0,         k_d;
      k_u,       k_hgs10,       -k_d,       k_r2;
        0,         k_u,          0,        -k_r2-k_d];

[V3, D3] = eig(K);
ev3 = diag(D3);


V3_0 = [1; 0; 0; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

A = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
B = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
C = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
D = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

V3_0 = [0; 1; 0; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

E = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
F = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
G = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
H = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

V3_0 = [0; 0; 1; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

I = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
J = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
K = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
L = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

V3_0 = [0; 0; 0; 1];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

M = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
N = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
O = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
P = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

%%
AA = @(T1,t2,T3,p) A(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
BB = @(T1,t2,T3,p) B(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
CC = @(T1,t2,T3,p) C(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
DD = @(T1,t2,T3,p) D(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
EE = @(T1,t2,T3,p) E(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
FF = @(T1,t2,T3,p) F(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
GG = @(T1,t2,T3,p) G(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
HH = @(T1,t2,T3,p) H(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
II = @(T1,t2,T3,p) I(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
JJ = @(T1,t2,T3,p) J(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
KK = @(T1,t2,T3,p) K(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
LL = @(T1,t2,T3,p) L(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
MM = @(T1,t2,T3,p) M(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
NN = @(T1,t2,T3,p) N(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
OO = @(T1,t2,T3,p) O(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
PP = @(T1,t2,T3,p) P(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);


%% ok let's try one inertial, one cone, and diffusion
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
dynParams1(1) = fitParamBnd('k_u',     0.0001, 0,  1, ''); %THESE CANT GO TO ZERO
dynParams1(2) = fitParamBnd('k_hgs10', 0.0001, 0,  1, ''); %THESE CANT GO TO ZERO
dynParams1(3) = fitParamBnd('k_r1',    0.0018, 0,  1, ''); % taken from T1 = 5.7 ps
dynParams1(4) = fitParamBnd('k_r2',    0.0050, 0,  1, ''); % taken from T1 = 2 psdyn(1) = aDArrayBnd(label,fun_array,ind_array,dynParams1,'fixed');
dyn(1) = aDArrayBnd(label,fun_array,ind_array,dynParams1,'fixed');

         
clear aRFoptions s;
aRFoptions.dt = 0.2;
aRFoptions.n_t = 96;
aRFoptions.n_zp = 2*aRFoptions.n_t;
aRFoptions.flag_rotating_frame = true;
aRFoptions.w_01_cm =     fitParamBnd('w_01_cm',   2059,   2020, 2080,'fixed');
aRFoptions.dw_sb_cm =    fitParamBnd('dw_sb_cm',   -11,  -15,  -10,  'fixed');
aRFoptions.mu01sq =      fitParamBnd('mu01sq',     0.08,   1e-5, 1.5, 'fixed');
aRFoptions.mu12sqRatio = fitParamBnd('mu12sqRatio',1.8,   0,    2.5, 'fixed');
aRFoptions.anh_cm =      fitParamBnd('anh_cm',      20,   10,   40,  'fixed');
aRFoptions.phase_deg =   fitParamBnd('phase_deg',    0,  -90,   90,  'fixed');
aRFoptions.dyn = dyn;
aRFoptions.dE_cm =       fitParamBnd('dE_cm',      700,   0,    1e5, 'fixed');
aRFoptions.temperature = fitParamBnd('temperature',292,   290,  505, 'fixed');
%aRFoptions.t2_array = [crop_data(scans).t2]./1000;
%aRFoptions.w1_in = crop_data(1).w1; %2320:0.5:2350;
%aRFoptions.w3_in = (crop_data(1).w3);%2296:2:2350;
aRFoptions.t2_array = [Aniso_19.t2]./1000;
aRFoptions.w1_in = w1; 
aRFoptions.w3_in = w3; %using Aniso data.w3 requires transpose
aRFoptions.useParallel=false;  
aRFoptions.nboot = 100;

% wobbling needs lsfOptions and aRFoptions 
% one angles and times are coming from a global fit of the anisotropy data
% at the center of the band
lsfParams(1) = fitParamBnd('Delta_cm',      20,     0.01,  50,     ''); %total linewidth
lsfParams(2) = fitParamBnd('theta0_deg',  24.7,     0,     180,    ''); %inertial cone   
lsfParams(3) = fitParamBnd('tr1',         4.42,     1,     100,    ''); %cone 1 time 
lsfParams(4) = fitParamBnd('theta1_deg',  69.4,     1,     179.999,''); %cone 1 angle
% lsfParams(5) = fitParamBnd('tr2',           70,     5,     200,    ''); %cone 2
% lsfParams(6) = fitParamBnd('theta2_deg',   30,   1,     180,    ''); %cone 2 angle
% lsfParams(7) = fitParamBnd('tr3',         100,   50,    500,    ''); %diffusive time or t3
% lsfParams(8) = fitParamBnd('theta3_deg',  90,   1,     180,    ''); %cone 3 angle

% best fit from bootstrapping CLS in clintonanalysis.m using
% RISDFitWobbling1inertial1cone1SSD.m . SSD terms come near the end as a
% CLS correction

aRFoptions.pol = 'para';
%aRFoptions.pol = 'perp';
aRFoptions.order = 1;

% the lsf needs params from the aRFoptions also to set up the numerical integration time
lsf = lsfRISDwobbling1inertial1coneNIBnd(lsfParams,'free',aRFoptions);

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
% s(1).dataMatrix = prepareGlobalFitData(crop_data(1,:))/(6e4);

%below should be for perp
aRFoptions.pol = 'perp';
lsf = lsfRISDwobbling1inertial1coneNIBnd(lsfParams,'free',aRFoptions); %have to reinitialize this so that the 'perp' get applied
fun_array = {oPerp}; 
dynOther = additionalDynamics(label,fun_array,ind_array);%has no free parameters (gets them from lsf)
aRFoptions.damping = lsf;
aRFoptions.dynOther = dynOther;
s(2) = aRFCO2Bnd(aRFoptions); %make a second s for perp spectra
% s(2).dataMatrix = prepareGlobalFitData(crop_data(2,:))/(1.5e4);

t=linspace(0,14,100);
p = lsf.copyParamValuesToParamStruct;
figure(31),plot(t,lsf.R.para(t,p)),title('FFCF R para and perp')
hold on
plot(t,lsf.R.perp(t,p))
hold off
figure(32),
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
%% chi by eye?
matrixPlot2DIR(m.s(1).simMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit',1); %sim
matrixPlot2DIR(m.s(2).simMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit',1); %sim

%% minimize error?
% m.globalFit %takes waaaaay too long to even start up the calculations?

%% output
% % 
% matrixPlot2DIR(m.s(1).simMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit'); %sim
% matrixPlot2DIR(m.s(2).simMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit'); %sim

% matrixPlot2DIR(m.s(1).dataMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit'); %real
% matrixPlot2DIR(m.s(2).dataMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit'); %real
%cool this thing works for once You just need to z-limit the spectra to
%see that they should have the same shape



%% time to get CLS values from the simMatrices?

zzzz_data = struct('w1',num2cell(1:28), 'w3',num2cell(1:28),'R',num2cell(1:28),'t2',num2cell(1:28));

for ii = 1 : 28
zzzz_data(ii).R = m.s(1).simMatrix(:,:,ii);
zzzz_data(ii).w1 = m.s(1).w1_in;
zzzz_data(ii).w3 = m.s(1).w3_in;
zzzz_data(ii).t2 = m.s(1).t2_array(ii).*1000;
end

d_zzzz = cropData(zzzz_data(1:end),w1,w3);
t2_array_zzzz =  [d_zzzz.t2]./1000;

zzxx_data = struct('w1',num2cell(1:28), 'w3',num2cell(1:28),'R',num2cell(1:28),'t2',num2cell(1:28));

for ii = 1 : 28
zzxx_data(ii).R = m.s(2).simMatrix(:,:,ii);
zzxx_data(ii).w1 = m.s(2).w1_in;
zzxx_data(ii).w3 = m.s(2).w3_in;
zzxx_data(ii).t2 = m.s(2).t2_array(ii).*1000;
end
  %% let's try and figure out how the ESA peak intensity changes as a function of time (normalized?)
  
%   for ii = 1 : length(t2_array_zzzz)
%       I_red_spara(ii) = zzzz_data(ii).R(6,22)./zzzz_data(1).R(6,22);
%       I_red_sperp(ii) = zzxx_data(ii).R(6,22)./zzxx_data(1).R(6,22);
%   
%       I_red_epara(ii) = crop_data(1,ii).R(6,21)./crop_data(1,1).R(6,21);
%       I_red_eperp(ii) = crop_data(2,ii).R(6,21)./crop_data(2,1).R(6,21);
%   end
%   
%   figure(501)
%   plot(t2_array_zzzz,I_red_spara,'ro')
%   hold on
%   plot(t2_array_zzzz,I_red_sperp,'rx')
%   plot(t2_array_zzzz,I_red_epara,'bo')
%   plot(t2_array_zzzz,I_red_eperp,'bx')
%       xlabel('t_2 (ps)')
%       ylabel('Normalized Intensity')
%       legend({'Sim para','Sim perp','Exp para','Exp perp'})
%       title('ESA Intensity')
%       xlim([0 210])
%       ylim([0 1.2])
%   hold off
  
  %% let's try and figure out how the GSB/SE peak intensity changes as a function of time (normalized?)
  
  for ii = 1 : length(t2_array_zzzz)
      I_blue_spara_norm(ii) = zzzz_data(ii).R(19,24)./zzzz_data(1).R(19,24);
      I_blue_sperp_norm(ii) = zzxx_data(ii).R(19,24)./zzxx_data(1).R(19,24);
  %center value becomes 70,70 for 1cm-1 resolution
%       I_blue_epara_norm(ii) = crop_data(1,ii).R(20,22)./crop_data(1,1).R(20,22);
%       I_blue_eperp_norm(ii) = crop_data(2,ii).R(20,22)./crop_data(2,1).R(20,22);
  end
  
  figure(301)
  plot(t2_array_zzzz,I_blue_spara_norm,'ro')
  hold on
  plot(t2_array_zzzz,I_blue_sperp_norm,'rx')
%   plot(t2_array_zzzz,I_blue_epara_norm,'bo')
%   plot(t2_array_zzzz,I_blue_eperp_norm,'bx')
      xlabel('t_2 (ps)')
      ylabel('Normalized Intensity')
      legend({'Sim para','Sim perp'})
      title('GSB/SE Intensity')
      xlim([0 12])
      ylim([0 1.2])
  hold off

  %% check the anisotropy and isotropic decays
% p=struct();

% p.k_u     = dynParams1(1).value;
% p.k_hgs10 = dynParams1(2).value;
% p.k_r     = dynParams1(3).value;
% p.dE_cm = 590;
% p.temperature = 290 ;
% p.Delta_cm   = lsfParams(1).value; %total linewidth
% % p.theta0_deg = lsfParams(2).value; %inertial cone   
% p.tr1        = lsfParams(3).value; %cone 1 time
% p.theta1_deg = lsfParams(4).value; %cone 1 angle
% p.tr2        = lsfParams(5).value; %cone 2
% p.theta2_deg = lsfParams(6).value; %cone 2 angle
% p.tr3        = lsfParams(7).value; %diffusive time or t3
% 
%   for ii = 1 : length(t2_array_zzzz)
%       I_blue_spara_unnorm(ii) = zzzz_data(ii).R(20,22);
%       I_blue_sperp_unnorm(ii) = zzxx_data(ii).R(20,22);
%       
%       I_blue_epara_unnorm(ii) = crop_data(1,ii).R(20,22);
%       I_blue_eperp_unnorm(ii) = crop_data(2,ii).R(20,22);
%   end
   
% 
% figure(2001)
% semilogx(t2_array_zzzz,(I_blue_epara_unnorm-I_blue_eperp_unnorm)./(I_blue_epara_unnorm+2.*I_blue_eperp_unnorm),'x')
% hold on
% semilogx(t2_array_zzzz,(I_blue_spara_unnorm-I_blue_sperp_unnorm)./(I_blue_spara_unnorm+2.*I_blue_sperp_unnorm),'o')
% semilogx(t2_array_zzzz,0.4.*L_l{2}(t2_array_zzzz,p))
%     xlabel('t_2 (ps)')
%     ylabel('anisotropy')
%     legend('Exp','Sim','L_2')
% hold off


% figure(2002)
% % semilogx(t2_array_zzzz,(I_blue_epara_norm+2.*I_blue_eperp_norm)./3,'rx')
% 
% semilogx(t2_array_zzzz,(lsf.R.para(t2_array_zzzz,p)+2.*lsf.R.perp(t2_array_zzzz,p))./3,'r')
% hold on
% semilogx(t2_array_zzzz,(I_blue_epara_norm+2.*I_blue_eperp_norm)./3,'ro')
% semilogx(t2_array_zzzz,(I_blue_spara_norm+2.*I_blue_sperp_norm)./3,'bo')
% semilogx(t2_array_zzzz,L_l{1}(t2_array_zzzz,p),'k')
%     xlabel('t_2 (ps)')
%     ylabel('isotropic decay')
% %     legend('Rpara+2Rperp','L_1')
% hold off


%%

% HWHM = 21;
% pfOptions.range1 = [2059 - HWHM   2059 + HWHM];
pfOptions.range1 = [2040 2082]; %range from paper1

pfOptions.range3 = [1990 2100]; 
view2DIRdata(zzzz_data(1),pfOptions.range1,pfOptions.range3);
pfOptions.fitfcn = fittype('- a1.*voigt(w,center,wg,wl) + a2.*voigt(w,center-anh,wg,wl) +c', ...
    'coeff',{'center','a1','a2','wg','wl','anh','c'},'indep',{'w'});
pfOptions.startpoint = [2059   1      1       7    7    15   0    ];
pfOptions.lb = [        2030   0      0       0    0    10   0    ];
pfOptions.ub = [        2100   Inf    Inf     100  100  30   1000 ];
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
% HWHM = 21;
% pfOptions.range1 = [2059 - HWHM   2059 + HWHM];
pfOptions.range1 = [2040 2082]; %range from paper1

pfOptions.range3 = [1990 2100];  
view2DIRdata(zzzz_data(1),pfOptions.range1,pfOptions.range3);
pfOptions.fitfcn = fittype('- a1.*voigt(w,center,wg,wl) + a2.*voigt(w,center-anh,wg,wl) + c', ...
    'coeff',{'center','a1','a2','wg','wl','anh','c'},'indep',{'w'});
pfOptions.startpoint = [2059   1      1       5    5    15   0    ];
pfOptions.lb = [        2030   0      0       0    0    10   0    ];
pfOptions.ub = [        2100   Inf    Inf     100  100  40   1000 ];
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
cfOptions(2).startpoint = [0.3 0.3 3    100 0.0 ];
cfOptions(2).lb =         [0   0   0.1  10  0    ];
cfOptions(2).ub =         [1   1   10   800 1    ];
cfOptions(2).fitfcn = fittype(@(a1,a2,t1,t2,c,x) ...
    a1.*exp(-x./t1) +  a2*exp(-x./t2) + c, ...
    'coeff', {'a1', 'a2', 't1','t2', 'c'}, 'indep', {'x'});
cfOptions(2).flag_plot = 0; %0=off 1=on 

cfOptions(3).startpoint = [0.2 0.2 0.2 3    10  100 ];
cfOptions(3).lb =         [0   0   0   0.1  5   50  ];
cfOptions(3).ub =         [1   1   1   5    50  Inf];
cfOptions(3).fitfcn = fittype(@(a1,a2,a3,t1,t2,t3,x) ...
      a1.*exp(-x./t1) +  a2*exp(-x./t2)+a3.*exp(-x./t3), ...
      'coeff', {'a1', 'a2','a3', 't1','t2','t3'}, 'indep', {'x'});
cfOptions(3).flag_plot = 0; %0=off 1=on

for ii = 1:3
    corr_struct_zzzz = corrFcnFit(t2_array_zzzz,c2_zzzz,c2_std_zzzz,cfOptions(ii));
    CFFit_zzzz(ii) = corr_struct_zzzz;
end
for ii = 1:3
    corr_struct_zzxx = corrFcnFit(t2_array_zzxx,c2_zzxx,c2_std_zzxx,cfOptions(ii));
    CFFit_zzxx(ii) = corr_struct_zzxx;
end

%%

figure(302) %biexponential
  set(gcf,'color','w')
%   plot(CFFit_zzzz(3).t2(:),(CFFit_zzzz(3).fitresult(CFFit_zzzz(3).t2(:)))); 
  hold on
  errorbar(CFFit_zzzz(2).t2(:),CFFit_zzzz(2).c2,CFFit_zzzz(2).c2_std,'bo'); 
  
%   plot(CFFit_zzxx(3).t2(:),(CFFit_zzxx(3).fitresult(CFFit_zzxx(3).t2(:)))); 
  errorbar(CFFit_zzxx(2).t2(:),CFFit_zzxx(2).c2,CFFit_zzxx(2).c2_std,'bx');
%   plot(t2_array_zzzz,lsf.R.para(t2_array_zzzz,p))
%   plot(t2_array_zzzz,lsf.R.perp(t2_array_zzzz,p))
plot(CFFit_CLS_LIFE006_para_19C.t2(:),CFFit_CLS_LIFE006_para_19C.fitresult(CFFit_CLS_LIFE006_para_19C.t2(:)),'k-');
%   errorbar(Bmimn2oRTv2avesepcal.analysis.zzzz(3).t2(:),Bmimn2oRTv2avesepcal.analysis.zzzz(3).c2(:),Bmimn2oRTv2avesepcal.analysis.zzzz(3).c2_std(:),'bo');
%   
  plot(CFFit_CLS_LIFE006_perp_19C.t2(:),CFFit_CLS_LIFE006_perp_19C.fitresult(CFFit_CLS_LIFE006_perp_19C.t2(:)),'k--');
%   errorbar(Bmimn2oRTv2avesepcal.analysis.zzxx(3).t2(:),Bmimn2oRTv2avesepcal.analysis.zzxx(3).c2(:),Bmimn2oRTv2avesepcal.analysis.zzxx(3).c2_std(:),'ro');
%  
%SSD on data points
% uh = (0.901.*exp(-(CFFit_zzzz(2).t2)./29.5)).*CFFit_zzzz(2).c2; %values for amp and tauSSD from bootstrapping
% uhhh = (0.901.*exp(-(CFFit_zzzz(2).t2)./29.5)).*CFFit_zzxx(2).c2; %CLS values in clintonanalysis.m

%SSD on R.para and Rperp
uh = (0.901.*exp(-(t2_array_zzzz)./29.5)).*lsf.R.para(t2_array_zzzz,p); %values for amp and tauSSD from bootstrapping
uhhh = (0.901.*exp(-(t2_array_zzzz)./29.5)).*lsf.R.perp(t2_array_zzzz,p); %CLS values in clintonanalysis.m

plot(CFFit_zzzz(3).t2,uh,'ro')
plot(CFFit_zzxx(3).t2,uhhh,'rx')
hold off
%% do some anisotropy

% for ii = 1:length(t2_array_zzzz)
% paramax(ii) = Bmimn2oRTv2avesepcal.g.para(ii).R(6,9);
% perpmax(ii) = Bmimn2oRTv2avesepcal.g.perp(ii).R(6,9); %all for real data
% 
% paracenter(ii) = Bmimn2oRTv2avesepcal.g.para(ii).R(5,7);
% perpcenter(ii) = Bmimn2oRTv2avesepcal.g.perp(ii).R(5,7);
% end
% 
% anisomax = (paramax - perpmax)./(paramax + 2.*perpmax);
% anisopeak = (paracenter - perpcenter)./(paracenter + 2.*perpcenter);
% t2 = t2_array_zzzz;

%% anisotropy of sim data
w1 = [2058 2061];
w3 = [2058 2061];
simpara = (prepareGlobalFitData(cropData(zzzz_data(1:end),w1,w3)));
simperp = (prepareGlobalFitData(cropData(zzxx_data(1:end),w1,w3)));

simaniso = (simpara - simperp)./(simpara + 2.*simperp);
simanisopeak = squeeze(simaniso(1,1,:))';
% figure(300)
%         [C,h]=contourf(cropData(zzzz_data(1),w1,w3).w1,cropData(zzzz_data(1),w1,w3).w3,simaniso(:,:,1),(0:0.02:0.7));
%         set(gcf,'PaperUnits','centimeters','PaperPosition',[5.0 5.0 7.5 5.0]);
%         title(sprintf('%0.1f ps',t2(1)));
%         xlabel('\omega_1/2\pic (cm^{-1})');
%         ylabel('\omega_3/2\pic (cm^{-1})');
%         clabel(C,h);

%%

figure(303)
plot([Aniso_19.t2]./1000,[Aniso_19.area],'x')
hold on
% plot(t2_array_zzzz,simanisopeak,'o')
plot(t2_array_zzzz,0.4.*lsf.L_l{2}(t2_array_zzzz,p),'k-')
xlim([0 12])
ylim([0 0.5])
% figure(300)
%         [C,h]=contourf(cropData(zzzz_data(1),w1,w3).w1,cropData(zzzz_data(1),w1,w3).w3,simaniso(:,:,1),(0:0.02:0.7));
%         set(gcf,'PaperUnits','centimeters','PaperPosition',[5.0 5.0 7.5 5.0]);
%         title(sprintf('%0.1f ps',t2(1)));
%         xlabel('\omega_1/2\pic (cm^{-1})');
%         ylabel('\omega_3/2\pic (cm^{-1})');
%         clabel(C,h);

%%






%% HERE FORWARD IS FOR Wobbling HL 1cone, 1SSD all in lineshape (g(t))
% close all
clear dyn

w1 = [1990:4:2130]; %change these!!!
w3 = [1990:2:2130];

% crop_data(1,1:5) = cropData(Bmimn2oRTv2avesepcal.zzzz([1 13 21 23 27]),w1,w3);
% crop_data(2,1:5) = cropData(Bmimn2oRTv2avesepcal.zzxx([1 13 21 23 27]),w1,w3);

syms k_u k_d k_hgs10 k_r1 k_r2 t


K = [-k_u,       k_r1,            k_d,        0;
        0,  -k_r1-k_hgs10-k_u,    0,         k_d;
      k_u,       k_hgs10,       -k_d,       k_r2;
        0,         k_u,          0,        -k_r2-k_d];

[V3, D3] = eig(K);
ev3 = diag(D3);


V3_0 = [1; 0; 0; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

A = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
B = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
C = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
D = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

V3_0 = [0; 1; 0; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

E = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
F = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
G = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
H = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

V3_0 = [0; 0; 1; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

I = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
J = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
K = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
L = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

V3_0 = [0; 0; 0; 1];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

M = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
N = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
O = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
P = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

%%
AA = @(T1,t2,T3,p) A(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
BB = @(T1,t2,T3,p) B(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
CC = @(T1,t2,T3,p) C(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
DD = @(T1,t2,T3,p) D(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
EE = @(T1,t2,T3,p) E(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
FF = @(T1,t2,T3,p) F(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
GG = @(T1,t2,T3,p) G(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
HH = @(T1,t2,T3,p) H(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
II = @(T1,t2,T3,p) I(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
JJ = @(T1,t2,T3,p) J(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
KK = @(T1,t2,T3,p) K(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
LL = @(T1,t2,T3,p) L(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
MM = @(T1,t2,T3,p) M(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
NN = @(T1,t2,T3,p) N(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
OO = @(T1,t2,T3,p) O(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
PP = @(T1,t2,T3,p) P(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);


%% ok let's try one inertial, one cone, and diffusion
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
dynParams1(1) = fitParamBnd('k_u',     1e-8,   0,Inf, ''); %THESE CANT GO TO ZERO
dynParams1(2) = fitParamBnd('k_hgs10', 1e-8,   0,  1, ''); %THESE CANT GO TO ZERO
dynParams1(3) = fitParamBnd('k_r1',    0.0018, 0,  1, ''); % taken from T1 = 5.7 ps
dynParams1(4) = fitParamBnd('k_r2',    0.0050, 0,  1, ''); % taken from T1 = 2 ps;
dyn(1) = aDArrayBnd(label,fun_array,ind_array,dynParams1,'fixed');

         
clear aRFoptions s;
aRFoptions.dt = 0.2;
aRFoptions.n_t = 96;
aRFoptions.n_zp = 2*aRFoptions.n_t;
aRFoptions.flag_rotating_frame = true;
aRFoptions.w_01_cm =     fitParamBnd('w_01_cm',   2059,  2020,  2100,'fixed'); %for low freq peak
aRFoptions.dw_sb_cm =    fitParamBnd('dw_sb_cm',   -11,   -19,   20,  'fixed');
aRFoptions.mu01sq =      fitParamBnd('mu01sq',    0.08,   1e-5, 1.5, 'fixed');
aRFoptions.mu12sqRatio = fitParamBnd('mu12sqRatio',1.8,   0,    2.5, 'fixed');
aRFoptions.anh_cm =      fitParamBnd('anh_cm',      20,   10,   40,  'fixed');
aRFoptions.phase_deg =   fitParamBnd('phase_deg',    0,  -90,   90,  'fixed');
aRFoptions.dyn = dyn;
aRFoptions.dE_cm =       fitParamBnd('dE_cm',      300,   0,    1e5, 'fixed'); %either 300 or 350
aRFoptions.temperature = fitParamBnd('temperature',292,   290,  505, 'fixed');
aRFoptions.w1_in = crop_data(1).w1; %2320:0.5:2350;
aRFoptions.w3_in = (crop_data(1).w3)';%2296:2:2350;
aRFoptions.t2_array = [para19.t2]./1000;
% aRFoptions.t2_array = [Aniso_19.t2]./1000;
% aRFoptions.w1_in = w1; 
% aRFoptions.w3_in = w3; %using Aniso data.w3 requires transpose
aRFoptions.useParallel=false;  
aRFoptions.nboot = 100;

% wobbling needs lsfOptions and aRFoptions 
% one angles and times are coming from a global fit of the anisotropy data
% at the center of the band
lsfParams(1) = fitParamBnd('Delta_cm',      18,     0.01,  50,     ''); %total linewidth 
% lsfParams(2) = fitParamBnd('theta0_deg',  24.7,     0,     180,    ''); %inertial cone   
lsfParams(3) = fitParamBnd('tr1',          2.6,     1,     100,    ''); %cone 1 time CAJ number
lsfParams(4) = fitParamBnd('theta1_deg',   108,     1,     179.999,''); %cone 1 angle CAJ number
lsfParams(5) = fitParamBnd('T2'        ,   2.8,     1,     50,     ''); %from FWHM to ps
lsfParams(6) = fitParamBnd('ampSSD'    ,  0.86,     0,     1 ,     ''); %from previous section of code (one cone)
lsfParams(7) = fitParamBnd('tauSSD'    ,    25,     1,    950,     ''); %from previous section of code (one cone)

% best fit from bootstrapping CLS in clintonanalysis.m using
% RISDFitWobbling1inertial1cone1SSD.m . SSD terms come near the end as a
% CLS correction

aRFoptions.pol = 'para';
%aRFoptions.pol = 'perp';
aRFoptions.order = 1;

% the lsf needs params from the aRFoptions also to set up the numerical integration time
lsf = lsfRISDwobblingHL1cone1SSDNIBnd(lsfParams,'free',aRFoptions);

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
s(1).dataMatrix = prepareGlobalFitData(crop_data(1,:))/(2e2);

%below should be for perp
aRFoptions.pol = 'perp';
lsf = lsfRISDwobblingHL1cone1SSDNIBnd(lsfParams,'free',aRFoptions); %have to reinitialize this so that the 'perp' get applied
fun_array = {oPerp}; 
dynOther = additionalDynamics(label,fun_array,ind_array);%has no free parameters (gets them from lsf)
aRFoptions.damping = lsf;
aRFoptions.dynOther = dynOther;
s(2) = aRFCO2Bnd(aRFoptions); %make a second s for perp spectra
s(2).dataMatrix = prepareGlobalFitData(crop_data(2,:))/(2e2);

t=linspace(0,15,1000);
p = lsf.copyParamValuesToParamStruct;
figure(31),plot(t,lsf.R.para(t,p)),title('FFCF R para and perp')
hold on
plot(t,lsf.R.perp(t,p))
hold off
figure(32),
for ii = 1:4
    subplot(4,1,ii)
    plot(t,lsf.L_l{ii}(t,p)),set(gca,'Ylim',[0,1]),title(['OCF L_{',num2str(ii),'}'])
end
% figure(500)
% hold on
% plot(t,0.4.*lsf.L_l{2}(t,p),'k')
% hold off
% plot([Aniso_19.t2]./1000,[Aniso_19.area],'x')

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
%% chi by eye?
matrixPlot2DIR(m.s(1).simMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit',1); %sim
matrixPlot2DIR(m.s(2).simMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit',1); %sim
% 
matrixPlot2DIR(m.s(1).dataMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit',1); %data
% matrixPlot2DIR(m.s(2).dataMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit',1); %data

%% minimize error?
% m.globalFit %takes waaaaay too long to even start up the calculations?

%% output
% % 
% matrixPlot2DIR(m.s(1).simMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit'); %sim
% matrixPlot2DIR(m.s(2).simMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit'); %sim

% matrixPlot2DIR(m.s(1).dataMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit'); %real
% matrixPlot2DIR(m.s(2).dataMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit'); %real
%cool this thing works for once You just need to z-limit the spectra to
%see that they should have the same shape



%% time to get CLS values from the simMatrices?

zzzz_data = struct('w1',num2cell(1:28), 'w3',num2cell(1:28),'R',num2cell(1:28),'t2',num2cell(1:28));

for ii = 1 : 28
zzzz_data(ii).R = m.s(1).simMatrix(:,:,ii);
zzzz_data(ii).w1 = m.s(1).w1_in;
zzzz_data(ii).w3 = m.s(1).w3_in;
zzzz_data(ii).t2 = m.s(1).t2_array(ii).*1000;
end

d_zzzz = cropData(zzzz_data(1:end),w1,w3);
t2_array_zzzz =  [d_zzzz.t2]./1000;

zzxx_data = struct('w1',num2cell(1:28), 'w3',num2cell(1:28),'R',num2cell(1:28),'t2',num2cell(1:28));

for ii = 1 : 28
zzxx_data(ii).R = m.s(2).simMatrix(:,:,ii);
zzxx_data(ii).w1 = m.s(2).w1_in;
zzxx_data(ii).w3 = m.s(2).w3_in;
zzxx_data(ii).t2 = m.s(2).t2_array(ii).*1000;
end
  %% let's try and figure out how the ESA peak intensity changes as a function of time (normalized?)
  
%   for ii = 1 : length(t2_array_zzzz)
%       I_red_spara(ii) = zzzz_data(ii).R(6,22)./zzzz_data(1).R(6,22);
%       I_red_sperp(ii) = zzxx_data(ii).R(6,22)./zzxx_data(1).R(6,22);
%   
%       I_red_epara(ii) = crop_data(1,ii).R(6,21)./crop_data(1,1).R(6,21);
%       I_red_eperp(ii) = crop_data(2,ii).R(6,21)./crop_data(2,1).R(6,21);
%   end
%   
%   figure(501)
%   plot(t2_array_zzzz,I_red_spara,'ro')
%   hold on
%   plot(t2_array_zzzz,I_red_sperp,'rx')
%   plot(t2_array_zzzz,I_red_epara,'bo')
%   plot(t2_array_zzzz,I_red_eperp,'bx')
%       xlabel('t_2 (ps)')
%       ylabel('Normalized Intensity')
%       legend({'Sim para','Sim perp','Exp para','Exp perp'})
%       title('ESA Intensity')
%       xlim([0 210])
%       ylim([0 1.2])
%   hold off
  
  %% let's try and figure out how the GSB/SE peak intensity changes as a function of time (normalized?)
  
  for ii = 1 : length(t2_array_zzzz)
      I_blue_spara_norm(ii) = zzzz_data(ii).R(30,18)./zzzz_data(1).R(30,18);
      I_blue_sperp_norm(ii) = zzxx_data(ii).R(30,18)./zzxx_data(1).R(30,18);

%       I_blue_epara_norm(ii) = crop_data(1,ii).R(20,22)./crop_data(1,1).R(20,22);
%       I_blue_eperp_norm(ii) = crop_data(2,ii).R(20,22)./crop_data(2,1).R(20,22);
  end
  
  figure(301)
  plot(t2_array_zzzz,I_blue_spara_norm,'ro')
  hold on
  plot(t2_array_zzzz,I_blue_sperp_norm,'rx')
%   plot(t2_array_zzzz,I_blue_epara_norm,'bo')
%   plot(t2_array_zzzz,I_blue_eperp_norm,'bx')
      xlabel('t_2 (ps)')
      ylabel('Normalized Intensity')
      legend({'Sim para','Sim perp'})
      title('GSB/SE Intensity')
      xlim([0 12])
      ylim([0 1.2])
  hold off

  %% check the anisotropy and isotropic decays
% p=struct();

% p.k_u     = dynParams1(1).value;
% p.k_hgs10 = dynParams1(2).value;
% p.k_r     = dynParams1(3).value;
% p.dE_cm = 590;
% p.temperature = 290 ;
% p.Delta_cm   = lsfParams(1).value; %total linewidth
% % p.theta0_deg = lsfParams(2).value; %inertial cone   
% p.tr1        = lsfParams(3).value; %cone 1 time
% p.theta1_deg = lsfParams(4).value; %cone 1 angle
% p.tr2        = lsfParams(5).value; %cone 2
% p.theta2_deg = lsfParams(6).value; %cone 2 angle
% p.tr3        = lsfParams(7).value; %diffusive time or t3
% 
%   for ii = 1 : length(t2_array_zzzz)
%       I_blue_spara_unnorm(ii) = zzzz_data(ii).R(20,22);
%       I_blue_sperp_unnorm(ii) = zzxx_data(ii).R(20,22);
%       
%       I_blue_epara_unnorm(ii) = crop_data(1,ii).R(20,22);
%       I_blue_eperp_unnorm(ii) = crop_data(2,ii).R(20,22);
%   end
   
% 
% figure(2001)
% semilogx(t2_array_zzzz,(I_blue_epara_unnorm-I_blue_eperp_unnorm)./(I_blue_epara_unnorm+2.*I_blue_eperp_unnorm),'x')
% hold on
% semilogx(t2_array_zzzz,(I_blue_spara_unnorm-I_blue_sperp_unnorm)./(I_blue_spara_unnorm+2.*I_blue_sperp_unnorm),'o')
% semilogx(t2_array_zzzz,0.4.*L_l{2}(t2_array_zzzz,p))
%     xlabel('t_2 (ps)')
%     ylabel('anisotropy')
%     legend('Exp','Sim','L_2')
% hold off


% figure(2002)
% % semilogx(t2_array_zzzz,(I_blue_epara_norm+2.*I_blue_eperp_norm)./3,'rx')
% 
% semilogx(t2_array_zzzz,(lsf.R.para(t2_array_zzzz,p)+2.*lsf.R.perp(t2_array_zzzz,p))./3,'r')
% hold on
% semilogx(t2_array_zzzz,(I_blue_epara_norm+2.*I_blue_eperp_norm)./3,'ro')
% semilogx(t2_array_zzzz,(I_blue_spara_norm+2.*I_blue_sperp_norm)./3,'bo')
% semilogx(t2_array_zzzz,L_l{1}(t2_array_zzzz,p),'k')
%     xlabel('t_2 (ps)')
%     ylabel('isotropic decay')
% %     legend('Rpara+2Rperp','L_1')
% hold off


%%

% HWHM = 21;
% pfOptions.range1 = [2059 - HWHM   2059 + HWHM];
pfOptions.range1 = [2040 2082]; %range from paper1

pfOptions.range3 = [1990 2100]; 
view2DIRdata(zzzz_data(1),pfOptions.range1,pfOptions.range3);
pfOptions.fitfcn = fittype('- a1.*voigt(w,center,wg,wl) + a2.*voigt(w,center-anh,wg,wl) +c', ...
    'coeff',{'center','a1','a2','wg','wl','anh','c'},'indep',{'w'});
pfOptions.startpoint = [2059   0.1    0.1     7    7    15   0    ];
pfOptions.lb = [        2030   0      0       0    0    10   0    ];
pfOptions.ub = [        2100   Inf    Inf     100  100  30   1000 ];
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
% HWHM = 21;
% pfOptions.range1 = [2059 - HWHM   2059 + HWHM];
pfOptions.range1 = [2040 2082]; %range from paper1

pfOptions.range3 = [1990 2100];  
view2DIRdata(zzzz_data(1),pfOptions.range1,pfOptions.range3);
pfOptions.fitfcn = fittype('- a1.*voigt(w,center,wg,wl) + a2.*voigt(w,center-anh,wg,wl) + c', ...
    'coeff',{'center','a1','a2','wg','wl','anh','c'},'indep',{'w'});
pfOptions.startpoint = [2059   0.1    0.1     5    5    15   0    ];
pfOptions.lb = [        2030   0      0       0    0    10   0    ];
pfOptions.ub = [        2100   Inf    Inf     100  100  40   1000 ];
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
% 
% cfOptions(1).startpoint = [0.3  3    0.0 ];
% cfOptions(1).lb =         [0    0.1  0    ];
% cfOptions(1).ub =         [1    500  1    ];
% cfOptions(1).fitfcn = fittype(@(a1,t1,c,x) ...
%     a1.*exp(-x./t1) + c, ...
%     'coeff', {'a1', 't1', 'c'}, 'indep', {'x'});
% cfOptions(1).flag_plot = 0; %0=off 1=on  
% 
cfOptions(2).startpoint = [0.3 0.3 3    100 0.0 ];
cfOptions(2).lb =         [0   0   0.1  10  0    ];
cfOptions(2).ub =         [1   1   10   800 1    ];
cfOptions(2).fitfcn = fittype(@(a1,a2,t1,t2,c,x) ...
    a1.*exp(-x./t1) +  a2*exp(-x./t2) + c, ...
    'coeff', {'a1', 'a2', 't1','t2', 'c'}, 'indep', {'x'});
cfOptions(2).flag_plot = 0; %0=off 1=on 

cfOptions(3).startpoint = [0.2 0.2 0.2 3    10  100 ];
cfOptions(3).lb =         [0   0   0   0.1  5   50  ];
cfOptions(3).ub =         [1   1   1   5    50  Inf];
cfOptions(3).fitfcn = fittype(@(a1,a2,a3,t1,t2,t3,x) ...
      a1.*exp(-x./t1) +  a2*exp(-x./t2)+a3.*exp(-x./t3), ...
      'coeff', {'a1', 'a2','a3', 't1','t2','t3'}, 'indep', {'x'});
cfOptions(3).flag_plot = 0; %0=off 1=on

for ii = 1:3
    corr_struct_zzzz = corrFcnFit(t2_array_zzzz,c2_zzzz,c2_std_zzzz,cfOptions(ii));
    CFFit_zzzz(ii) = corr_struct_zzzz;
end
for ii = 1:3
    corr_struct_zzxx = corrFcnFit(t2_array_zzxx,c2_zzxx,c2_std_zzxx,cfOptions(ii));
    CFFit_zzxx(ii) = corr_struct_zzxx;
end

%%

figure(302) %biexponential
  set(gcf,'color','w')
%   plot(CFFit_zzzz(3).t2(:),(CFFit_zzzz(3).fitresult(CFFit_zzzz(3).t2(:)))); 
  hold on
  errorbar(CFFit_zzzz(2).t2(:),CFFit_zzzz(2).c2,CFFit_zzzz(2).c2_std,'bo'); 
  
%   plot(CFFit_zzxx(3).t2(:),(CFFit_zzxx(3).fitresult(CFFit_zzxx(3).t2(:)))); 
  errorbar(CFFit_zzxx(2).t2(:),CFFit_zzxx(2).c2,CFFit_zzxx(2).c2_std,'bx');
%   plot(t2_array_zzzz,lsf.R.para(t2_array_zzzz,p))
%   plot(t2_array_zzzz,lsf.R.perp(t2_array_zzzz,p))
plot(CFFit_CLS_LIFE006_para_19C.t2(:),CFFit_CLS_LIFE006_para_19C.fitresult(CFFit_CLS_LIFE006_para_19C.t2(:)),'k-');
%   errorbar(Bmimn2oRTv2avesepcal.analysis.zzzz(3).t2(:),Bmimn2oRTv2avesepcal.analysis.zzzz(3).c2(:),Bmimn2oRTv2avesepcal.analysis.zzzz(3).c2_std(:),'bo');
%   
  plot(CFFit_CLS_LIFE006_perp_19C.t2(:),CFFit_CLS_LIFE006_perp_19C.fitresult(CFFit_CLS_LIFE006_perp_19C.t2(:)),'k--');
%   errorbar(Bmimn2oRTv2avesepcal.analysis.zzxx(3).t2(:),Bmimn2oRTv2avesepcal.analysis.zzxx(3).c2(:),Bmimn2oRTv2avesepcal.analysis.zzxx(3).c2_std(:),'ro');
%  
%SSD on data points
% uh = (0.901.*exp(-(CFFit_zzzz(2).t2)./29.5)).*CFFit_zzzz(2).c2; %values for amp and tauSSD from bootstrapping
% uhhh = (0.901.*exp(-(CFFit_zzzz(2).t2)./29.5)).*CFFit_zzxx(2).c2; %CLS values in clintonanalysis.m

%SSD on R.para and Rperp
% uh = (0.901.*exp(-(t2_array_zzzz)./29.5)).*lsf.R.para(t2_array_zzzz,p); %values for amp and tauSSD from bootstrapping
% uhhh = (0.901.*exp(-(t2_array_zzzz)./29.5)).*lsf.R.perp(t2_array_zzzz,p); %CLS values in clintonanalysis.m

% plot(CFFit_zzzz(3).t2,uh,'ro')
% plot(CFFit_zzxx(3).t2,uhhh,'rx')
xlabel('t_2 (ps)')
ylabel('c_2')
hold off
%% do some anisotropy

% for ii = 1:length(t2_array_zzzz)
% paramax(ii) = Bmimn2oRTv2avesepcal.g.para(ii).R(6,9);
% perpmax(ii) = Bmimn2oRTv2avesepcal.g.perp(ii).R(6,9); %all for real data
% 
% paracenter(ii) = Bmimn2oRTv2avesepcal.g.para(ii).R(5,7);
% perpcenter(ii) = Bmimn2oRTv2avesepcal.g.perp(ii).R(5,7);
% end
% 
% anisomax = (paramax - perpmax)./(paramax + 2.*perpmax);
% anisopeak = (paracenter - perpcenter)./(paracenter + 2.*perpcenter);
% t2 = t2_array_zzzz;

%% anisotropy of sim data
w1a = [2058 2062];
w3a = [2058 2062];
simpara = (prepareGlobalFitData(cropData(zzzz_data(1:end),w1a,w3a)));
simperp = (prepareGlobalFitData(cropData(zzxx_data(1:end),w1a,w3a)));

simaniso = (simpara - simperp)./(simpara + 2.*simperp);
simanisopeak = squeeze(simaniso(1,1,:))';
% figure(300)
%         [C,h]=contourf(cropData(zzzz_data(1),w1,w3).w1,cropData(zzzz_data(1),w1,w3).w3,simaniso(:,:,1),(0:0.02:0.7));
%         set(gcf,'PaperUnits','centimeters','PaperPosition',[5.0 5.0 7.5 5.0]);
%         title(sprintf('%0.1f ps',t2(1)));
%         xlabel('\omega_1/2\pic (cm^{-1})');
%         ylabel('\omega_3/2\pic (cm^{-1})');
%         clabel(C,h);

%%

figure(303);clf
plot([Aniso_19.t2]./1000,[Aniso_19.area],'kx')
hold on
plot(t2_array_zzzz,simanisopeak,'bo')
plot(t2_array_zzzz,0.4.*lsf.L_l{2}(t2_array_zzzz,p),'b-')
xlim([0 12])
ylim([0 0.5])
xlabel('t_2 (ps)')
ylabel('Anisotropy')
box off
hold off
% figure(300)
%         [C,h]=contourf(cropData(zzzz_data(1),w1,w3).w1,cropData(zzzz_data(1),w1,w3).w3,simaniso(:,:,1),(0:0.02:0.7));
%         set(gcf,'PaperUnits','centimeters','PaperPosition',[5.0 5.0 7.5 5.0]);
%         title(sprintf('%0.1f ps',t2(1)));
%         xlabel('\omega_1/2\pic (cm^{-1})');
%         ylabel('\omega_3/2\pic (cm^{-1})');
%         clabel(C,h);
%%






%% HERE FORWARD IS FOR Wobbling HL 1cone, 1SSD all in lineshape (g(t)) using anisotropy fit values
% close all
clear dyn

w1 = [1990:4:2130]; %change these!!!
w3 = [1990:2:2130];

% crop_data(1,1:5) = cropData(Bmimn2oRTv2avesepcal.zzzz([1 13 21 23 27]),w1,w3);
% crop_data(2,1:5) = cropData(Bmimn2oRTv2avesepcal.zzxx([1 13 21 23 27]),w1,w3);

syms k_u k_d k_hgs10 k_r1 k_r2 t


K = [-k_u,       k_r1,            k_d,        0;
        0,  -k_r1-k_hgs10-k_u,    0,         k_d;
      k_u,       k_hgs10,       -k_d,       k_r2;
        0,         k_u,          0,        -k_r2-k_d];

[V3, D3] = eig(K);
ev3 = diag(D3);


V3_0 = [1; 0; 0; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

A = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
B = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
C = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
D = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

V3_0 = [0; 1; 0; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

E = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
F = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
G = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
H = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

V3_0 = [0; 0; 1; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

I = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
J = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
K = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
L = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

V3_0 = [0; 0; 0; 1];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

M = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
N = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
O = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
P = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

%%
AA = @(T1,t2,T3,p) A(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
BB = @(T1,t2,T3,p) B(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
CC = @(T1,t2,T3,p) C(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
DD = @(T1,t2,T3,p) D(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
EE = @(T1,t2,T3,p) E(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
FF = @(T1,t2,T3,p) F(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
GG = @(T1,t2,T3,p) G(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
HH = @(T1,t2,T3,p) H(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
II = @(T1,t2,T3,p) I(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
JJ = @(T1,t2,T3,p) J(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
KK = @(T1,t2,T3,p) K(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
LL = @(T1,t2,T3,p) L(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
MM = @(T1,t2,T3,p) M(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
NN = @(T1,t2,T3,p) N(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
OO = @(T1,t2,T3,p) O(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
PP = @(T1,t2,T3,p) P(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);


%% ok let's try one inertial, one cone, and diffusion
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
dynParams1(1) = fitParamBnd('k_u',     1e-8,   0,Inf, ''); %THESE CANT GO TO ZERO
dynParams1(2) = fitParamBnd('k_hgs10', 1e-8,   0,  1, ''); %THESE CANT GO TO ZERO
dynParams1(3) = fitParamBnd('k_r1',    0.0018, 0,  1, ''); % taken from T1 = 5.7 ps
dynParams1(4) = fitParamBnd('k_r2',    0.0050, 0,  1, ''); % taken from T1 = 2 ps;
dyn(1) = aDArrayBnd(label,fun_array,ind_array,dynParams1,'fixed');

         
clear aRFoptions s;
aRFoptions.dt = 0.2;
aRFoptions.n_t = 96;
aRFoptions.n_zp = 2*aRFoptions.n_t;
aRFoptions.flag_rotating_frame = true;
aRFoptions.w_01_cm =     fitParamBnd('w_01_cm',   2059,  2020,  2100,'fixed'); %for low freq peak
aRFoptions.dw_sb_cm =    fitParamBnd('dw_sb_cm',   -11,   -19,   20,  'fixed');
aRFoptions.mu01sq =      fitParamBnd('mu01sq',    0.08,   1e-5, 1.5, 'fixed');
aRFoptions.mu12sqRatio = fitParamBnd('mu12sqRatio',1.8,   0,    2.5, 'fixed');
aRFoptions.anh_cm =      fitParamBnd('anh_cm',      20,   10,   40,  'fixed');
aRFoptions.phase_deg =   fitParamBnd('phase_deg',    0,  -90,   90,  'fixed');
aRFoptions.dyn = dyn;
aRFoptions.dE_cm =       fitParamBnd('dE_cm',      400,   0,    1e5, 'fixed'); %either 300 or 350
aRFoptions.temperature = fitParamBnd('temperature',292,   290,  505, 'fixed');
aRFoptions.w1_in = crop_data(1).w1; %2320:0.5:2350;
aRFoptions.w3_in = (crop_data(1).w3)';%2296:2:2350;
aRFoptions.t2_array = [para19.t2]./1000;
% aRFoptions.t2_array = [Aniso_19.t2]./1000;
% aRFoptions.w1_in = w1; 
% aRFoptions.w3_in = w3; %using Aniso data.w3 requires transpose
aRFoptions.useParallel=false;  
aRFoptions.nboot = 100;

% wobbling needs lsfOptions and aRFoptions 
% one angles and times are coming from a global fit of the anisotropy data
% at the center of the band
lsfParams(1) = fitParamBnd('Delta_cm',      17,     0.01,  50,     ''); %total linewidth 
% lsfParams(2) = fitParamBnd('theta0_deg',  24.7,     0,     180,    ''); %inertial cone   
lsfParams(3) = fitParamBnd('tr1',         20.0,     1,     100,    ''); %cone 1 time CAJ number
lsfParams(4) = fitParamBnd('theta1_deg',    70,     1,     179.999,''); %cone 1 angle CAJ number
lsfParams(5) = fitParamBnd('T2'        ,  1.20,     1,     50,     ''); %from averaging T1 times
lsfParams(6) = fitParamBnd('ampSSD'    ,  0.90,     0,     1 ,     ''); %from previous section of code (one cone)
lsfParams(7) = fitParamBnd('tauSSD'    ,    15,     1,    950,     ''); %from previous section of code (one cone)

% best fit from bootstrapping CLS in clintonanalysis.m using
% RISDFitWobbling1inertial1cone1SSD.m . SSD terms come near the end as a
% CLS correction

aRFoptions.pol = 'para';
%aRFoptions.pol = 'perp';
aRFoptions.order = 1;

% the lsf needs params from the aRFoptions also to set up the numerical integration time
lsf = lsfRISDwobblingHL1cone1SSDNIBnd(lsfParams,'free',aRFoptions);

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
s(1).dataMatrix = prepareGlobalFitData(crop_data(1,:))/(2e2);

%below should be for perp
aRFoptions.pol = 'perp';
lsf = lsfRISDwobblingHL1cone1SSDNIBnd(lsfParams,'free',aRFoptions); %have to reinitialize this so that the 'perp' get applied
fun_array = {oPerp}; 
dynOther = additionalDynamics(label,fun_array,ind_array);%has no free parameters (gets them from lsf)
aRFoptions.damping = lsf;
aRFoptions.dynOther = dynOther;
s(2) = aRFCO2Bnd(aRFoptions); %make a second s for perp spectra
s(2).dataMatrix = prepareGlobalFitData(crop_data(2,:))/(2e2);

t=linspace(0,15,1000);
p = lsf.copyParamValuesToParamStruct;
figure(31),plot(t,lsf.R.para(t,p)),title('FFCF R para and perp')
hold on
plot(t,lsf.R.perp(t,p))
hold off
figure(32),
for ii = 1:4
    subplot(4,1,ii)
    plot(t,lsf.L_l{ii}(t,p)),set(gca,'Ylim',[0,1]),title(['OCF L_{',num2str(ii),'}'])
end
% figure(500)
% hold on
% plot(t,0.4.*lsf.L_l{2}(t,p),'k')
% hold off
% plot([Aniso_19.t2]./1000,[Aniso_19.area],'x')

m = multiMeasurement(s);

%% simulate via multiMeasurement

%simulated
tic
m.s(1) = m.s(1).calcSpectrum(m.s(1).p0);
% matrixPlot2DIR(s(1).simMatrix,s(1).w1_in,s(1).w3_in,s(1).t2_array,[5 6],'height_px',100,'zlimit',1) %sim

%perpendicular
m.s(2) = m.s(2).calcSpectrum(m.s(2).p0);
% matrixPlot2DIR(s(2).simMatrix,s(2).w1_in,s(2).w3_in,s(2).t2_array,[5 6],'height_px',100,'zlimit',1) %sim
toc
%% chi by eye?
matrixPlot2DIR(m.s(1).simMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit',1); %sim
matrixPlot2DIR(m.s(2).simMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit',1); %sim
% 
matrixPlot2DIR(m.s(1).dataMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit',1); %data
% matrixPlot2DIR(m.s(2).dataMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit',1); %data

%% minimize error?
% m.globalFit %takes waaaaay too long to even start up the calculations?

%% output
% % 
% matrixPlot2DIR(m.s(1).simMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit'); %sim
% matrixPlot2DIR(m.s(2).simMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit'); %sim

% matrixPlot2DIR(m.s(1).dataMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit'); %real
% matrixPlot2DIR(m.s(2).dataMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit'); %real
%cool this thing works for once You just need to z-limit the spectra to
%see that they should have the same shape



%% time to get CLS values from the simMatrices?

zzzz_data = struct('w1',num2cell(1:28), 'w3',num2cell(1:28),'R',num2cell(1:28),'t2',num2cell(1:28));

for ii = 1 : 28
zzzz_data(ii).R = m.s(1).simMatrix(:,:,ii);
zzzz_data(ii).w1 = m.s(1).w1_in;
zzzz_data(ii).w3 = m.s(1).w3_in;
zzzz_data(ii).t2 = m.s(1).t2_array(ii).*1000;
end

d_zzzz = cropData(zzzz_data(1:end),w1,w3);
t2_array_zzzz =  [d_zzzz.t2]./1000;

zzxx_data = struct('w1',num2cell(1:28), 'w3',num2cell(1:28),'R',num2cell(1:28),'t2',num2cell(1:28));

for ii = 1 : 28
zzxx_data(ii).R = m.s(2).simMatrix(:,:,ii);
zzxx_data(ii).w1 = m.s(2).w1_in;
zzxx_data(ii).w3 = m.s(2).w3_in;
zzxx_data(ii).t2 = m.s(2).t2_array(ii).*1000;
end
  %% let's try and figure out how the ESA peak intensity changes as a function of time (normalized?)
  
%   for ii = 1 : length(t2_array_zzzz)
%       I_red_spara(ii) = zzzz_data(ii).R(6,22)./zzzz_data(1).R(6,22);
%       I_red_sperp(ii) = zzxx_data(ii).R(6,22)./zzxx_data(1).R(6,22);
%   
%       I_red_epara(ii) = crop_data(1,ii).R(6,21)./crop_data(1,1).R(6,21);
%       I_red_eperp(ii) = crop_data(2,ii).R(6,21)./crop_data(2,1).R(6,21);
%   end
%   
%   figure(501)
%   plot(t2_array_zzzz,I_red_spara,'ro')
%   hold on
%   plot(t2_array_zzzz,I_red_sperp,'rx')
%   plot(t2_array_zzzz,I_red_epara,'bo')
%   plot(t2_array_zzzz,I_red_eperp,'bx')
%       xlabel('t_2 (ps)')
%       ylabel('Normalized Intensity')
%       legend({'Sim para','Sim perp','Exp para','Exp perp'})
%       title('ESA Intensity')
%       xlim([0 210])
%       ylim([0 1.2])
%   hold off
  
  %% let's try and figure out how the GSB/SE peak intensity changes as a function of time (normalized?)
  
  for ii = 1 : length(t2_array_zzzz)
      I_blue_spara_norm(ii) = zzzz_data(ii).R(30,18)./zzzz_data(1).R(30,18);
      I_blue_sperp_norm(ii) = zzxx_data(ii).R(30,18)./zzxx_data(1).R(30,18);

%       I_blue_epara_norm(ii) = crop_data(1,ii).R(20,22)./crop_data(1,1).R(20,22);
%       I_blue_eperp_norm(ii) = crop_data(2,ii).R(20,22)./crop_data(2,1).R(20,22);
  end
  
  figure(301)
  plot(t2_array_zzzz,I_blue_spara_norm,'ro')
  hold on
  plot(t2_array_zzzz,I_blue_sperp_norm,'rx')
%   plot(t2_array_zzzz,I_blue_epara_norm,'bo')
%   plot(t2_array_zzzz,I_blue_eperp_norm,'bx')
      xlabel('t_2 (ps)')
      ylabel('Normalized Intensity')
      legend({'Sim para','Sim perp'})
      title('GSB/SE Intensity')
      xlim([0 12])
      ylim([0 1.2])
  hold off

  %% check the anisotropy and isotropic decays
% p=struct();

% p.k_u     = dynParams1(1).value;
% p.k_hgs10 = dynParams1(2).value;
% p.k_r     = dynParams1(3).value;
% p.dE_cm = 590;
% p.temperature = 290 ;
% p.Delta_cm   = lsfParams(1).value; %total linewidth
% % p.theta0_deg = lsfParams(2).value; %inertial cone   
% p.tr1        = lsfParams(3).value; %cone 1 time
% p.theta1_deg = lsfParams(4).value; %cone 1 angle
% p.tr2        = lsfParams(5).value; %cone 2
% p.theta2_deg = lsfParams(6).value; %cone 2 angle
% p.tr3        = lsfParams(7).value; %diffusive time or t3
% 
%   for ii = 1 : length(t2_array_zzzz)
%       I_blue_spara_unnorm(ii) = zzzz_data(ii).R(20,22);
%       I_blue_sperp_unnorm(ii) = zzxx_data(ii).R(20,22);
%       
%       I_blue_epara_unnorm(ii) = crop_data(1,ii).R(20,22);
%       I_blue_eperp_unnorm(ii) = crop_data(2,ii).R(20,22);
%   end
   
% 
% figure(2001)
% semilogx(t2_array_zzzz,(I_blue_epara_unnorm-I_blue_eperp_unnorm)./(I_blue_epara_unnorm+2.*I_blue_eperp_unnorm),'x')
% hold on
% semilogx(t2_array_zzzz,(I_blue_spara_unnorm-I_blue_sperp_unnorm)./(I_blue_spara_unnorm+2.*I_blue_sperp_unnorm),'o')
% semilogx(t2_array_zzzz,0.4.*L_l{2}(t2_array_zzzz,p))
%     xlabel('t_2 (ps)')
%     ylabel('anisotropy')
%     legend('Exp','Sim','L_2')
% hold off


% figure(2002)
% % semilogx(t2_array_zzzz,(I_blue_epara_norm+2.*I_blue_eperp_norm)./3,'rx')
% 
% semilogx(t2_array_zzzz,(lsf.R.para(t2_array_zzzz,p)+2.*lsf.R.perp(t2_array_zzzz,p))./3,'r')
% hold on
% semilogx(t2_array_zzzz,(I_blue_epara_norm+2.*I_blue_eperp_norm)./3,'ro')
% semilogx(t2_array_zzzz,(I_blue_spara_norm+2.*I_blue_sperp_norm)./3,'bo')
% semilogx(t2_array_zzzz,L_l{1}(t2_array_zzzz,p),'k')
%     xlabel('t_2 (ps)')
%     ylabel('isotropic decay')
% %     legend('Rpara+2Rperp','L_1')
% hold off


%%

% HWHM = 21;
% pfOptions.range1 = [2059 - HWHM   2059 + HWHM];
pfOptions.range1 = [2040 2082]; %range from paper1

pfOptions.range3 = [1990 2100]; 
view2DIRdata(zzzz_data(1),pfOptions.range1,pfOptions.range3);
pfOptions.fitfcn = fittype('- a1.*voigt(w,center,wg,wl) + a2.*voigt(w,center-anh,wg,wl) +c', ...
    'coeff',{'center','a1','a2','wg','wl','anh','c'},'indep',{'w'});
pfOptions.startpoint = [2059   0.1    0.1     7    7    15   0    ];
pfOptions.lb = [        2030   0      0       0    0    10   0    ];
pfOptions.ub = [        2100   Inf    Inf     100  100  30   1000 ];
pfOptions.zp_factor = 2;
pfOptions.flag_plot = 0;
pfOptions.estimateArea = 1;

d_zzzz = cropData(zzzz_data(1:end),pfOptions.range1,pfOptions.range3);
t2_array_zzzz =  [d_zzzz.t2]./1000;

%%
peakFit_zzzz = extractMaxima(d_zzzz,pfOptions);

%%
maxMatrix_zzzz = CLSMaxMatrix(peakFit_zzzz);
[CLS_zzzz,c2_zzzz,c2_std_zzzz] = fitCLS(d_zzzz,maxMatrix_zzzz,1);

%% ZZXX fitting
% HWHM = 21;
% pfOptions.range1 = [2059 - HWHM   2059 + HWHM];
pfOptions.range1 = [2040 2082]; %range from paper1

pfOptions.range3 = [1990 2100];  
view2DIRdata(zzzz_data(1),pfOptions.range1,pfOptions.range3);
pfOptions.fitfcn = fittype('- a1.*voigt(w,center,wg,wl) + a2.*voigt(w,center-anh,wg,wl) + c', ...
    'coeff',{'center','a1','a2','wg','wl','anh','c'},'indep',{'w'});
pfOptions.startpoint = [2059   0.1    0.1     5    5    15   0    ];
pfOptions.lb = [        2030   0      0       0    0    10   0    ];
pfOptions.ub = [        2100   Inf    Inf     100  100  40   1000 ];
pfOptions.zp_factor = 2;
pfOptions.flag_plot = 0;
pfOptions.estimateArea = 1;

d_zzxx = cropData(zzxx_data(1:end),pfOptions.range1,pfOptions.range3);
t2_array_zzxx =  [d_zzxx.t2]./1000;

%%
peakFit_zzxx = extractMaxima(d_zzxx,pfOptions);
%%
maxMatrix_zzxx = CLSMaxMatrix(peakFit_zzxx);
[CLS_zzxx,c2_zzxx,c2_std_zzxx] = fitCLS(d_zzxx,maxMatrix_zzxx,1);

%% Fitting Corr.Funct
% close all
% 
% cfOptions(1).startpoint = [0.3  3    0.0 ];
% cfOptions(1).lb =         [0    0.1  0    ];
% cfOptions(1).ub =         [1    500  1    ];
% cfOptions(1).fitfcn = fittype(@(a1,t1,c,x) ...
%     a1.*exp(-x./t1) + c, ...
%     'coeff', {'a1', 't1', 'c'}, 'indep', {'x'});
% cfOptions(1).flag_plot = 0; %0=off 1=on  
% 
cfOptions(2).startpoint = [0.3 0.3 3    100 0.0 ];
cfOptions(2).lb =         [0   0   0.1  10  0    ];
cfOptions(2).ub =         [1   1   10   800 1    ];
cfOptions(2).fitfcn = fittype(@(a1,a2,t1,t2,c,x) ...
    a1.*exp(-x./t1) +  a2*exp(-x./t2) + c, ...
    'coeff', {'a1', 'a2', 't1','t2', 'c'}, 'indep', {'x'});
cfOptions(2).flag_plot = 0; %0=off 1=on 

cfOptions(3).startpoint = [0.2 0.2 0.2 3    10  100 ];
cfOptions(3).lb =         [0   0   0   0.1  5   50  ];
cfOptions(3).ub =         [1   1   1   5    50  Inf];
cfOptions(3).fitfcn = fittype(@(a1,a2,a3,t1,t2,t3,x) ...
      a1.*exp(-x./t1) +  a2*exp(-x./t2)+a3.*exp(-x./t3), ...
      'coeff', {'a1', 'a2','a3', 't1','t2','t3'}, 'indep', {'x'});
cfOptions(3).flag_plot = 0; %0=off 1=on

for ii = 2:3
    corr_struct_zzzz = corrFcnFit(t2_array_zzzz,c2_zzzz,c2_std_zzzz,cfOptions(ii));
    CFFit_zzzz(ii) = corr_struct_zzzz;
end
for ii = 2:3
    corr_struct_zzxx = corrFcnFit(t2_array_zzxx,c2_zzxx,c2_std_zzxx,cfOptions(ii));
    CFFit_zzxx(ii) = corr_struct_zzxx;
end

%%

figure(302) %biexponential
  set(gcf,'color','w')
%   plot(CFFit_zzzz(3).t2(:),(CFFit_zzzz(3).fitresult(CFFit_zzzz(3).t2(:)))); 
  hold on
  errorbar(CFFit_zzzz(2).t2(:),CFFit_zzzz(2).c2,CFFit_zzzz(2).c2_std,'bo'); 
  
%   plot(CFFit_zzxx(3).t2(:),(CFFit_zzxx(3).fitresult(CFFit_zzxx(3).t2(:)))); 
  errorbar(CFFit_zzxx(2).t2(:),CFFit_zzxx(2).c2,CFFit_zzxx(2).c2_std,'bx');
%   plot(t2_array_zzzz,lsf.R.para(t2_array_zzzz,p))
%   plot(t2_array_zzzz,lsf.R.perp(t2_array_zzzz,p))
plot(CFFit_CLS_LIFE006_para_19C.t2(:),CFFit_CLS_LIFE006_para_19C.fitresult(CFFit_CLS_LIFE006_para_19C.t2(:)),'k-');
%   errorbar(Bmimn2oRTv2avesepcal.analysis.zzzz(3).t2(:),Bmimn2oRTv2avesepcal.analysis.zzzz(3).c2(:),Bmimn2oRTv2avesepcal.analysis.zzzz(3).c2_std(:),'bo');
%   
  plot(CFFit_CLS_LIFE006_perp_19C.t2(:),CFFit_CLS_LIFE006_perp_19C.fitresult(CFFit_CLS_LIFE006_perp_19C.t2(:)),'k--');
%   errorbar(Bmimn2oRTv2avesepcal.analysis.zzxx(3).t2(:),Bmimn2oRTv2avesepcal.analysis.zzxx(3).c2(:),Bmimn2oRTv2avesepcal.analysis.zzxx(3).c2_std(:),'ro');
%  
%SSD on data points
% uh = (0.901.*exp(-(CFFit_zzzz(2).t2)./29.5)).*CFFit_zzzz(2).c2; %values for amp and tauSSD from bootstrapping
% uhhh = (0.901.*exp(-(CFFit_zzzz(2).t2)./29.5)).*CFFit_zzxx(2).c2; %CLS values in clintonanalysis.m

%SSD on R.para and Rperp
% uh = (0.901.*exp(-(t2_array_zzzz)./29.5)).*lsf.R.para(t2_array_zzzz,p); %values for amp and tauSSD from bootstrapping
% uhhh = (0.901.*exp(-(t2_array_zzzz)./29.5)).*lsf.R.perp(t2_array_zzzz,p); %CLS values in clintonanalysis.m

% plot(CFFit_zzzz(3).t2,uh,'ro')
% plot(CFFit_zzxx(3).t2,uhhh,'rx')
xlabel('t_2 (ps)')
ylabel('c_2')
hold off
%% do some anisotropy

% for ii = 1:length(t2_array_zzzz)
% paramax(ii) = Bmimn2oRTv2avesepcal.g.para(ii).R(6,9);
% perpmax(ii) = Bmimn2oRTv2avesepcal.g.perp(ii).R(6,9); %all for real data
% 
% paracenter(ii) = Bmimn2oRTv2avesepcal.g.para(ii).R(5,7);
% perpcenter(ii) = Bmimn2oRTv2avesepcal.g.perp(ii).R(5,7);
% end
% 
% anisomax = (paramax - perpmax)./(paramax + 2.*perpmax);
% anisopeak = (paracenter - perpcenter)./(paracenter + 2.*perpcenter);
% t2 = t2_array_zzzz;

%% anisotropy of sim data
w1a = [2058 2062];
w3a = [2058 2062];
simpara = (prepareGlobalFitData(cropData(zzzz_data(1:end),w1a,w3a)));
simperp = (prepareGlobalFitData(cropData(zzxx_data(1:end),w1a,w3a)));

simaniso = (simpara - simperp)./(simpara + 2.*simperp);
simanisopeak = squeeze(simaniso(1,1,:))';
% figure(300)
%         [C,h]=contourf(cropData(zzzz_data(1),w1,w3).w1,cropData(zzzz_data(1),w1,w3).w3,simaniso(:,:,1),(0:0.02:0.7));
%         set(gcf,'PaperUnits','centimeters','PaperPosition',[5.0 5.0 7.5 5.0]);
%         title(sprintf('%0.1f ps',t2(1)));
%         xlabel('\omega_1/2\pic (cm^{-1})');
%         ylabel('\omega_3/2\pic (cm^{-1})');
%         clabel(C,h);

%%

figure(303);clf
plot([Aniso_19.t2]./1000,[Aniso_19.area],'kx')
hold on
plot(t2_array_zzzz,simanisopeak,'bo')
plot(t2_array_zzzz,0.4.*lsf.L_l{2}(t2_array_zzzz,p),'b-')
xlim([0 12])
ylim([0 0.5])
xlabel('t_2 (ps)')
ylabel('Anisotropy')
box off
hold off
% figure(300)
%         [C,h]=contourf(cropData(zzzz_data(1),w1,w3).w1,cropData(zzzz_data(1),w1,w3).w3,simaniso(:,:,1),(0:0.02:0.7));
%         set(gcf,'PaperUnits','centimeters','PaperPosition',[5.0 5.0 7.5 5.0]);
%         title(sprintf('%0.1f ps',t2(1)));
%         xlabel('\omega_1/2\pic (cm^{-1})');
%         ylabel('\omega_3/2\pic (cm^{-1})');
%         clabel(C,h);


%% HERE FORWARD IS FOR Wobbling HL 1cone, 1SSD all in lineshape (g(t)) using anisotropy fit values 
% high freq peak is the center!
% close all
clear dyn

w1 = [1990:4:2130]; %change these!!!
w3 = [1990:2:2130];

% crop_data(1,1:5) = cropData(Bmimn2oRTv2avesepcal.zzzz([1 13 21 23 27]),w1,w3);
% crop_data(2,1:5) = cropData(Bmimn2oRTv2avesepcal.zzxx([1 13 21 23 27]),w1,w3);

syms k_u k_d k_hgs10 k_r1 k_r2 t


K = [-k_u,       k_r1,            k_d,        0;
        0,  -k_r1-k_hgs10-k_u,    0,         k_d;
      k_u,       k_hgs10,       -k_d,       k_r2;
        0,         k_u,          0,        -k_r2-k_d];

[V3, D3] = eig(K);
ev3 = diag(D3);


V3_0 = [1; 0; 0; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

A = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
B = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
C = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
D = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

V3_0 = [0; 1; 0; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

E = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
F = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
G = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
H = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

V3_0 = [0; 0; 1; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

I = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
J = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
K = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
L = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

V3_0 = [0; 0; 0; 1];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

M = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
N = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
O = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
P = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

%%
AA = @(T1,t2,T3,p) A(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
BB = @(T1,t2,T3,p) B(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
CC = @(T1,t2,T3,p) C(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
DD = @(T1,t2,T3,p) D(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
EE = @(T1,t2,T3,p) E(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
FF = @(T1,t2,T3,p) F(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
GG = @(T1,t2,T3,p) G(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
HH = @(T1,t2,T3,p) H(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
II = @(T1,t2,T3,p) I(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
JJ = @(T1,t2,T3,p) J(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
KK = @(T1,t2,T3,p) K(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
LL = @(T1,t2,T3,p) L(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
MM = @(T1,t2,T3,p) M(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
NN = @(T1,t2,T3,p) N(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
OO = @(T1,t2,T3,p) O(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
PP = @(T1,t2,T3,p) P(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);


%% ok let's try one inertial, one cone, and diffusion
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
dynParams1(1) = fitParamBnd('k_u',     0.01,   0,Inf, ''); %THESE CANT GO TO ZERO
dynParams1(2) = fitParamBnd('k_hgs10', 1e-8,   0,  1, ''); %THESE CANT GO TO ZERO
dynParams1(3) = fitParamBnd('k_r1',    0.0050, 0,  1, ''); % taken from T1 = 5.7 ps
dynParams1(4) = fitParamBnd('k_r2',    0.0018, 0,  1, ''); % taken from T1 = 2 ps;
dyn(1) = aDArrayBnd(label,fun_array,ind_array,dynParams1,'fixed');

         
clear aRFoptions s;
aRFoptions.dt = 0.2;
aRFoptions.n_t = 96;
aRFoptions.n_zp = 2*aRFoptions.n_t;
aRFoptions.flag_rotating_frame = true;
aRFoptions.w_01_cm =     fitParamBnd('w_01_cm',   2059,  2020,  2100,'fixed'); %for low freq peak
aRFoptions.dw_sb_cm =    fitParamBnd('dw_sb_cm',   -11,   -19,   20,  'fixed');
aRFoptions.mu01sq =      fitParamBnd('mu01sq',    0.08,   1e-5, 1.5, 'fixed');
aRFoptions.mu12sqRatio = fitParamBnd('mu12sqRatio',1.8,   0,    2.5, 'fixed');
aRFoptions.anh_cm =      fitParamBnd('anh_cm',      20,   10,   40,  'fixed');
aRFoptions.phase_deg =   fitParamBnd('phase_deg',    0,  -90,   90,  'fixed');
aRFoptions.dyn = dyn;
aRFoptions.dE_cm =       fitParamBnd('dE_cm',      390,   0,    1e5, 'fixed'); %either 300 or 350
aRFoptions.temperature = fitParamBnd('temperature',292,   290,  505, 'fixed');
aRFoptions.w1_in = crop_data(1).w1; %2320:0.5:2350;
aRFoptions.w3_in = (crop_data(1).w3)';%2296:2:2350;
aRFoptions.t2_array = [para19.t2]./1000;
% aRFoptions.t2_array = [Aniso_19.t2]./1000;
% aRFoptions.w1_in = w1; 
% aRFoptions.w3_in = w3; %using Aniso data.w3 requires transpose
aRFoptions.useParallel=false;  
aRFoptions.nboot = 100;

% wobbling needs lsfOptions and aRFoptions 
% one angles and times are coming from a global fit of the anisotropy data
% at the center of the band
lsfParams(1) = fitParamBnd('Delta_cm',      17,     0.01,  50,     ''); %total linewidth 
% lsfParams(2) = fitParamBnd('theta0_deg',  24.7,     0,     180,    ''); %inertial cone   
lsfParams(3) = fitParamBnd('tr1',         20.0,     1,     100,    ''); %cone 1 time CAJ number
lsfParams(4) = fitParamBnd('theta1_deg',    70,     1,     179.999,''); %cone 1 angle CAJ number
lsfParams(5) = fitParamBnd('T2'        ,  1.20,     1,     50,     ''); %from averaging T1 times
lsfParams(6) = fitParamBnd('ampSSD'    ,  0.90,     0,     1 ,     ''); %from previous section of code (one cone)
lsfParams(7) = fitParamBnd('tauSSD'    ,    15,     1,    950,     ''); %from previous section of code (one cone)

% best fit from bootstrapping CLS in clintonanalysis.m using
% RISDFitWobbling1inertial1cone1SSD.m . SSD terms come near the end as a
% CLS correction

aRFoptions.pol = 'para';
%aRFoptions.pol = 'perp';
aRFoptions.order = 1;

% the lsf needs params from the aRFoptions also to set up the numerical integration time
lsf = lsfRISDwobblingHL1cone1SSDNIBnd(lsfParams,'free',aRFoptions);

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
s(1).dataMatrix = prepareGlobalFitData(crop_data(1,:))/(2e2);

%below should be for perp
aRFoptions.pol = 'perp';
lsf = lsfRISDwobblingHL1cone1SSDNIBnd(lsfParams,'free',aRFoptions); %have to reinitialize this so that the 'perp' get applied
fun_array = {oPerp}; 
dynOther = additionalDynamics(label,fun_array,ind_array);%has no free parameters (gets them from lsf)
aRFoptions.damping = lsf;
aRFoptions.dynOther = dynOther;
s(2) = aRFCO2Bnd(aRFoptions); %make a second s for perp spectra
s(2).dataMatrix = prepareGlobalFitData(crop_data(2,:))/(2e2);

t=linspace(0,15,1000);
p = lsf.copyParamValuesToParamStruct;
figure(31),plot(t,lsf.R.para(t,p)),title('FFCF R para and perp')
hold on
plot(t,lsf.R.perp(t,p))
hold off
figure(32),
for ii = 1:4
    subplot(4,1,ii)
    plot(t,lsf.L_l{ii}(t,p)),set(gca,'Ylim',[0,1]),title(['OCF L_{',num2str(ii),'}'])
end
% figure(500)
% hold on
% plot(t,0.4.*lsf.L_l{2}(t,p),'k')
% hold off
% plot([Aniso_19.t2]./1000,[Aniso_19.area],'x')

m = multiMeasurement(s);

%% simulate via multiMeasurement

%simulated
tic
m.s(1) = m.s(1).calcSpectrum(m.s(1).p0);
% matrixPlot2DIR(s(1).simMatrix,s(1).w1_in,s(1).w3_in,s(1).t2_array,[5 6],'height_px',100,'zlimit',1) %sim

%perpendicular
m.s(2) = m.s(2).calcSpectrum(m.s(2).p0);
% matrixPlot2DIR(s(2).simMatrix,s(2).w1_in,s(2).w3_in,s(2).t2_array,[5 6],'height_px',100,'zlimit',1) %sim
toc
%% chi by eye?
matrixPlot2DIR(m.s(1).simMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit',1); %sim
matrixPlot2DIR(m.s(2).simMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit',1); %sim
% 
matrixPlot2DIR(m.s(1).dataMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit',1); %data
% matrixPlot2DIR(m.s(2).dataMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit',1); %data

%% minimize error?
% m.globalFit %takes waaaaay too long to even start up the calculations?

%% output
% % 
% matrixPlot2DIR(m.s(1).simMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit'); %sim
% matrixPlot2DIR(m.s(2).simMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit'); %sim

% matrixPlot2DIR(m.s(1).dataMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit'); %real
% matrixPlot2DIR(m.s(2).dataMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit'); %real
%cool this thing works for once You just need to z-limit the spectra to
%see that they should have the same shape



%% time to get CLS values from the simMatrices?

zzzz_data = struct('w1',num2cell(1:28), 'w3',num2cell(1:28),'R',num2cell(1:28),'t2',num2cell(1:28));

for ii = 1 : 28
zzzz_data(ii).R = m.s(1).simMatrix(:,:,ii);
zzzz_data(ii).w1 = m.s(1).w1_in;
zzzz_data(ii).w3 = m.s(1).w3_in;
zzzz_data(ii).t2 = m.s(1).t2_array(ii).*1000;
end

d_zzzz = cropData(zzzz_data(1:end),w1,w3);
t2_array_zzzz =  [d_zzzz.t2]./1000;

zzxx_data = struct('w1',num2cell(1:28), 'w3',num2cell(1:28),'R',num2cell(1:28),'t2',num2cell(1:28));

for ii = 1 : 28
zzxx_data(ii).R = m.s(2).simMatrix(:,:,ii);
zzxx_data(ii).w1 = m.s(2).w1_in;
zzxx_data(ii).w3 = m.s(2).w3_in;
zzxx_data(ii).t2 = m.s(2).t2_array(ii).*1000;
end
  %% let's try and figure out how the ESA peak intensity changes as a function of time (normalized?)
  
%   for ii = 1 : length(t2_array_zzzz)
%       I_red_spara(ii) = zzzz_data(ii).R(6,22)./zzzz_data(1).R(6,22);
%       I_red_sperp(ii) = zzxx_data(ii).R(6,22)./zzxx_data(1).R(6,22);
%   
%       I_red_epara(ii) = crop_data(1,ii).R(6,21)./crop_data(1,1).R(6,21);
%       I_red_eperp(ii) = crop_data(2,ii).R(6,21)./crop_data(2,1).R(6,21);
%   end
%   
%   figure(501)
%   plot(t2_array_zzzz,I_red_spara,'ro')
%   hold on
%   plot(t2_array_zzzz,I_red_sperp,'rx')
%   plot(t2_array_zzzz,I_red_epara,'bo')
%   plot(t2_array_zzzz,I_red_eperp,'bx')
%       xlabel('t_2 (ps)')
%       ylabel('Normalized Intensity')
%       legend({'Sim para','Sim perp','Exp para','Exp perp'})
%       title('ESA Intensity')
%       xlim([0 210])
%       ylim([0 1.2])
%   hold off
  
  %% let's try and figure out how the GSB/SE peak intensity changes as a function of time (normalized?)
  
  for ii = 1 : length(t2_array_zzzz)
      I_blue_spara_norm(ii) = zzzz_data(ii).R(30,18)./zzzz_data(1).R(30,18);
      I_blue_sperp_norm(ii) = zzxx_data(ii).R(30,18)./zzxx_data(1).R(30,18);

%       I_blue_epara_norm(ii) = crop_data(1,ii).R(20,22)./crop_data(1,1).R(20,22);
%       I_blue_eperp_norm(ii) = crop_data(2,ii).R(20,22)./crop_data(2,1).R(20,22);
  end
  
  figure(301)
  plot(t2_array_zzzz,I_blue_spara_norm,'ro')
  hold on
  plot(t2_array_zzzz,I_blue_sperp_norm,'rx')
%   plot(t2_array_zzzz,I_blue_epara_norm,'bo')
%   plot(t2_array_zzzz,I_blue_eperp_norm,'bx')
      xlabel('t_2 (ps)')
      ylabel('Normalized Intensity')
      legend({'Sim para','Sim perp'})
      title('GSB/SE Intensity')
      xlim([0 12])
      ylim([0 1.2])
  hold off

  %% check the anisotropy and isotropic decays
% p=struct();

% p.k_u     = dynParams1(1).value;
% p.k_hgs10 = dynParams1(2).value;
% p.k_r     = dynParams1(3).value;
% p.dE_cm = 590;
% p.temperature = 290 ;
% p.Delta_cm   = lsfParams(1).value; %total linewidth
% % p.theta0_deg = lsfParams(2).value; %inertial cone   
% p.tr1        = lsfParams(3).value; %cone 1 time
% p.theta1_deg = lsfParams(4).value; %cone 1 angle
% p.tr2        = lsfParams(5).value; %cone 2
% p.theta2_deg = lsfParams(6).value; %cone 2 angle
% p.tr3        = lsfParams(7).value; %diffusive time or t3
% 
%   for ii = 1 : length(t2_array_zzzz)
%       I_blue_spara_unnorm(ii) = zzzz_data(ii).R(20,22);
%       I_blue_sperp_unnorm(ii) = zzxx_data(ii).R(20,22);
%       
%       I_blue_epara_unnorm(ii) = crop_data(1,ii).R(20,22);
%       I_blue_eperp_unnorm(ii) = crop_data(2,ii).R(20,22);
%   end
   
% 
% figure(2001)
% semilogx(t2_array_zzzz,(I_blue_epara_unnorm-I_blue_eperp_unnorm)./(I_blue_epara_unnorm+2.*I_blue_eperp_unnorm),'x')
% hold on
% semilogx(t2_array_zzzz,(I_blue_spara_unnorm-I_blue_sperp_unnorm)./(I_blue_spara_unnorm+2.*I_blue_sperp_unnorm),'o')
% semilogx(t2_array_zzzz,0.4.*L_l{2}(t2_array_zzzz,p))
%     xlabel('t_2 (ps)')
%     ylabel('anisotropy')
%     legend('Exp','Sim','L_2')
% hold off


% figure(2002)
% % semilogx(t2_array_zzzz,(I_blue_epara_norm+2.*I_blue_eperp_norm)./3,'rx')
% 
% semilogx(t2_array_zzzz,(lsf.R.para(t2_array_zzzz,p)+2.*lsf.R.perp(t2_array_zzzz,p))./3,'r')
% hold on
% semilogx(t2_array_zzzz,(I_blue_epara_norm+2.*I_blue_eperp_norm)./3,'ro')
% semilogx(t2_array_zzzz,(I_blue_spara_norm+2.*I_blue_sperp_norm)./3,'bo')
% semilogx(t2_array_zzzz,L_l{1}(t2_array_zzzz,p),'k')
%     xlabel('t_2 (ps)')
%     ylabel('isotropic decay')
% %     legend('Rpara+2Rperp','L_1')
% hold off


%%

% HWHM = 21;
% pfOptions.range1 = [2059 - HWHM   2059 + HWHM];
pfOptions.range1 = [2040 2082]; %range from paper1

pfOptions.range3 = [1990 2100]; 
view2DIRdata(zzzz_data(1),pfOptions.range1,pfOptions.range3);
pfOptions.fitfcn = fittype('- a1.*voigt(w,center,wg,wl) + a2.*voigt(w,center-anh,wg,wl) +c', ...
    'coeff',{'center','a1','a2','wg','wl','anh','c'},'indep',{'w'});
pfOptions.startpoint = [2059   0.1    0.1     7    7    15   0    ];
pfOptions.lb = [        2030   0      0       0    0    10   0    ];
pfOptions.ub = [        2100   Inf    Inf     100  100  30   1000 ];
pfOptions.zp_factor = 2;
pfOptions.flag_plot = 0;
pfOptions.estimateArea = 1;

d_zzzz = cropData(zzzz_data(1:end),pfOptions.range1,pfOptions.range3);
t2_array_zzzz =  [d_zzzz.t2]./1000;

%%
peakFit_zzzz = extractMaxima(d_zzzz,pfOptions);

%%
maxMatrix_zzzz = CLSMaxMatrix(peakFit_zzzz);
[CLS_zzzz,c2_zzzz,c2_std_zzzz] = fitCLS(d_zzzz,maxMatrix_zzzz,1);

%% ZZXX fitting
% HWHM = 21;
% pfOptions.range1 = [2059 - HWHM   2059 + HWHM];
pfOptions.range1 = [2040 2082]; %range from paper1

pfOptions.range3 = [1990 2100];  
view2DIRdata(zzzz_data(1),pfOptions.range1,pfOptions.range3);
pfOptions.fitfcn = fittype('- a1.*voigt(w,center,wg,wl) + a2.*voigt(w,center-anh,wg,wl) + c', ...
    'coeff',{'center','a1','a2','wg','wl','anh','c'},'indep',{'w'});
pfOptions.startpoint = [2059   0.1    0.1     5    5    15   0    ];
pfOptions.lb = [        2030   0      0       0    0    10   0    ];
pfOptions.ub = [        2100   Inf    Inf     100  100  40   1000 ];
pfOptions.zp_factor = 2;
pfOptions.flag_plot = 0;
pfOptions.estimateArea = 1;

d_zzxx = cropData(zzxx_data(1:end),pfOptions.range1,pfOptions.range3);
t2_array_zzxx =  [d_zzxx.t2]./1000;

%%
peakFit_zzxx = extractMaxima(d_zzxx,pfOptions);
%%
maxMatrix_zzxx = CLSMaxMatrix(peakFit_zzxx);
[CLS_zzxx,c2_zzxx,c2_std_zzxx] = fitCLS(d_zzxx,maxMatrix_zzxx,1);

%% Fitting Corr.Funct
% close all
% 
% cfOptions(1).startpoint = [0.3  3    0.0 ];
% cfOptions(1).lb =         [0    0.1  0    ];
% cfOptions(1).ub =         [1    500  1    ];
% cfOptions(1).fitfcn = fittype(@(a1,t1,c,x) ...
%     a1.*exp(-x./t1) + c, ...
%     'coeff', {'a1', 't1', 'c'}, 'indep', {'x'});
% cfOptions(1).flag_plot = 0; %0=off 1=on  
% 
cfOptions(2).startpoint = [0.3 0.3 3    100 0.0 ];
cfOptions(2).lb =         [0   0   0.1  10  0    ];
cfOptions(2).ub =         [1   1   10   800 1    ];
cfOptions(2).fitfcn = fittype(@(a1,a2,t1,t2,c,x) ...
    a1.*exp(-x./t1) +  a2*exp(-x./t2) + c, ...
    'coeff', {'a1', 'a2', 't1','t2', 'c'}, 'indep', {'x'});
cfOptions(2).flag_plot = 0; %0=off 1=on 

cfOptions(3).startpoint = [0.2 0.2 0.2 3    10  100 ];
cfOptions(3).lb =         [0   0   0   0.1  5   50  ];
cfOptions(3).ub =         [1   1   1   5    50  Inf];
cfOptions(3).fitfcn = fittype(@(a1,a2,a3,t1,t2,t3,x) ...
      a1.*exp(-x./t1) +  a2*exp(-x./t2)+a3.*exp(-x./t3), ...
      'coeff', {'a1', 'a2','a3', 't1','t2','t3'}, 'indep', {'x'});
cfOptions(3).flag_plot = 0; %0=off 1=on

for ii = 2:3
    corr_struct_zzzz = corrFcnFit(t2_array_zzzz,c2_zzzz,c2_std_zzzz,cfOptions(ii));
    CFFit_zzzz(ii) = corr_struct_zzzz;
end
for ii = 2:3
    corr_struct_zzxx = corrFcnFit(t2_array_zzxx,c2_zzxx,c2_std_zzxx,cfOptions(ii));
    CFFit_zzxx(ii) = corr_struct_zzxx;
end

%%

figure(302) %biexponential
  set(gcf,'color','w')
%   plot(CFFit_zzzz(3).t2(:),(CFFit_zzzz(3).fitresult(CFFit_zzzz(3).t2(:)))); 
  hold on
  errorbar(CFFit_zzzz(2).t2(:),CFFit_zzzz(2).c2,CFFit_zzzz(2).c2_std,'bo'); 
  
%   plot(CFFit_zzxx(3).t2(:),(CFFit_zzxx(3).fitresult(CFFit_zzxx(3).t2(:)))); 
  errorbar(CFFit_zzxx(2).t2(:),CFFit_zzxx(2).c2,CFFit_zzxx(2).c2_std,'bx');
%   plot(t2_array_zzzz,lsf.R.para(t2_array_zzzz,p))
%   plot(t2_array_zzzz,lsf.R.perp(t2_array_zzzz,p))
plot(CFFit_CLS_LIFE006_para_19C.t2(:),CFFit_CLS_LIFE006_para_19C.fitresult(CFFit_CLS_LIFE006_para_19C.t2(:)),'k-');
%   errorbar(Bmimn2oRTv2avesepcal.analysis.zzzz(3).t2(:),Bmimn2oRTv2avesepcal.analysis.zzzz(3).c2(:),Bmimn2oRTv2avesepcal.analysis.zzzz(3).c2_std(:),'bo');
%   
  plot(CFFit_CLS_LIFE006_perp_19C.t2(:),CFFit_CLS_LIFE006_perp_19C.fitresult(CFFit_CLS_LIFE006_perp_19C.t2(:)),'k--');
%   errorbar(Bmimn2oRTv2avesepcal.analysis.zzxx(3).t2(:),Bmimn2oRTv2avesepcal.analysis.zzxx(3).c2(:),Bmimn2oRTv2avesepcal.analysis.zzxx(3).c2_std(:),'ro');
%  
%SSD on data points
% uh = (0.901.*exp(-(CFFit_zzzz(2).t2)./29.5)).*CFFit_zzzz(2).c2; %values for amp and tauSSD from bootstrapping
% uhhh = (0.901.*exp(-(CFFit_zzzz(2).t2)./29.5)).*CFFit_zzxx(2).c2; %CLS values in clintonanalysis.m

%SSD on R.para and Rperp
% uh = (0.901.*exp(-(t2_array_zzzz)./29.5)).*lsf.R.para(t2_array_zzzz,p); %values for amp and tauSSD from bootstrapping
% uhhh = (0.901.*exp(-(t2_array_zzzz)./29.5)).*lsf.R.perp(t2_array_zzzz,p); %CLS values in clintonanalysis.m

% plot(CFFit_zzzz(3).t2,uh,'ro')
% plot(CFFit_zzxx(3).t2,uhhh,'rx')
xlabel('t_2 (ps)')
ylabel('c_2')
hold off
%% do some anisotropy

% for ii = 1:length(t2_array_zzzz)
% paramax(ii) = Bmimn2oRTv2avesepcal.g.para(ii).R(6,9);
% perpmax(ii) = Bmimn2oRTv2avesepcal.g.perp(ii).R(6,9); %all for real data
% 
% paracenter(ii) = Bmimn2oRTv2avesepcal.g.para(ii).R(5,7);
% perpcenter(ii) = Bmimn2oRTv2avesepcal.g.perp(ii).R(5,7);
% end
% 
% anisomax = (paramax - perpmax)./(paramax + 2.*perpmax);
% anisopeak = (paracenter - perpcenter)./(paracenter + 2.*perpcenter);
% t2 = t2_array_zzzz;

%% anisotropy of sim data
w1a = [2058 2062];
w3a = [2058 2062];
simpara = (prepareGlobalFitData(cropData(zzzz_data(1:end),w1a,w3a)));
simperp = (prepareGlobalFitData(cropData(zzxx_data(1:end),w1a,w3a)));

simaniso = (simpara - simperp)./(simpara + 2.*simperp);
simanisopeak = squeeze(simaniso(1,1,:))';
% figure(300)
%         [C,h]=contourf(cropData(zzzz_data(1),w1,w3).w1,cropData(zzzz_data(1),w1,w3).w3,simaniso(:,:,1),(0:0.02:0.7));
%         set(gcf,'PaperUnits','centimeters','PaperPosition',[5.0 5.0 7.5 5.0]);
%         title(sprintf('%0.1f ps',t2(1)));
%         xlabel('\omega_1/2\pic (cm^{-1})');
%         ylabel('\omega_3/2\pic (cm^{-1})');
%         clabel(C,h);

%%

figure(303);clf
plot([Aniso_19.t2]./1000,[Aniso_19.area],'kx')
hold on
plot(t2_array_zzzz,simanisopeak,'bo')
plot(t2_array_zzzz,0.4.*lsf.L_l{2}(t2_array_zzzz,p),'b-')
xlim([0 12])
ylim([0 0.5])
xlabel('t_2 (ps)')
ylabel('Anisotropy')
box off
hold off
% figure(300)
%         [C,h]=contourf(cropData(zzzz_data(1),w1,w3).w1,cropData(zzzz_data(1),w1,w3).w3,simaniso(:,:,1),(0:0.02:0.7));
%         set(gcf,'PaperUnits','centimeters','PaperPosition',[5.0 5.0 7.5 5.0]);
%         title(sprintf('%0.1f ps',t2(1)));
%         xlabel('\omega_1/2\pic (cm^{-1})');
%         ylabel('\omega_3/2\pic (cm^{-1})');
%         clabel(C,h);



%% Just going to do a diffusive orientation cone (from bootstrapping)

%% HERE FORWARD IS FOR Wobbling HL 1cone, 1SSD all in lineshape (g(t)) using anisotropy fit values 
% high freq peak is the center!
% close all
clear dyn

w1 = [1990:4:2130]; %change these!!!
w3 = [1990:2:2130];

% crop_data(1,1:5) = cropData(Bmimn2oRTv2avesepcal.zzzz([1 13 21 23 27]),w1,w3);
% crop_data(2,1:5) = cropData(Bmimn2oRTv2avesepcal.zzxx([1 13 21 23 27]),w1,w3);

syms k_u k_d k_hgs10 k_r1 k_r2 t


K = [-k_u,       k_r1,            k_d,        0;
        0,  -k_r1-k_hgs10-k_u,    0,         k_d;
      k_u,       k_hgs10,       -k_d,       k_r2;
        0,         k_u,          0,        -k_r2-k_d];

[V3, D3] = eig(K);
ev3 = diag(D3);


V3_0 = [1; 0; 0; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

A = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
B = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
C = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
D = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

V3_0 = [0; 1; 0; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

E = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
F = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
G = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
H = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

V3_0 = [0; 0; 1; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

I = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
J = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
K = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
L = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

V3_0 = [0; 0; 0; 1];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

M = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
N = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
O = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);
P = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_r1,k_r2,t]);

%%
AA = @(T1,t2,T3,p) A(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
BB = @(T1,t2,T3,p) B(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
CC = @(T1,t2,T3,p) C(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
DD = @(T1,t2,T3,p) D(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
EE = @(T1,t2,T3,p) E(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
FF = @(T1,t2,T3,p) F(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
GG = @(T1,t2,T3,p) G(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
HH = @(T1,t2,T3,p) H(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
II = @(T1,t2,T3,p) I(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
JJ = @(T1,t2,T3,p) J(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
KK = @(T1,t2,T3,p) K(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
LL = @(T1,t2,T3,p) L(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
MM = @(T1,t2,T3,p) M(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
NN = @(T1,t2,T3,p) N(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
OO = @(T1,t2,T3,p) O(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);
PP = @(T1,t2,T3,p) P(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_r1,p.k_r2,t2);


%% ok let's try one inertial, one cone, and diffusion
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
dynParams1(1) = fitParamBnd('k_u',     0.01,   0,Inf, ''); %THESE CANT GO TO ZERO
dynParams1(2) = fitParamBnd('k_hgs10', 1e-8,   0,  1, ''); %THESE CANT GO TO ZERO
dynParams1(3) = fitParamBnd('k_r1',    0.0050, 0,  1, ''); % taken from T1 = 5.7 ps
dynParams1(4) = fitParamBnd('k_r2',    0.0018, 0,  1, ''); % taken from T1 = 2 ps;
dyn(1) = aDArrayBnd(label,fun_array,ind_array,dynParams1,'fixed');

         
clear aRFoptions s;
aRFoptions.dt = 0.2;
aRFoptions.n_t = 96;
aRFoptions.n_zp = 2*aRFoptions.n_t;
aRFoptions.flag_rotating_frame = true;
aRFoptions.w_01_cm =     fitParamBnd('w_01_cm',   2059,  2020,  2100,'fixed'); %for low freq peak
aRFoptions.dw_sb_cm =    fitParamBnd('dw_sb_cm',   -11,   -19,   20,  'fixed');
aRFoptions.mu01sq =      fitParamBnd('mu01sq',    0.08,   1e-5, 1.5, 'fixed');
aRFoptions.mu12sqRatio = fitParamBnd('mu12sqRatio',1.8,   0,    2.5, 'fixed');
aRFoptions.anh_cm =      fitParamBnd('anh_cm',      20,   10,   40,  'fixed');
aRFoptions.phase_deg =   fitParamBnd('phase_deg',    0,  -90,   90,  'fixed');
aRFoptions.dyn = dyn;
aRFoptions.dE_cm =       fitParamBnd('dE_cm',      390,   0,    1e5, 'fixed'); %either 300 or 350
aRFoptions.temperature = fitParamBnd('temperature',292,   290,  505, 'fixed');
aRFoptions.w1_in = crop_data(1).w1; %2320:0.5:2350;
aRFoptions.w3_in = (crop_data(1).w3)';%2296:2:2350;
aRFoptions.t2_array = [para19.t2]./1000;
% aRFoptions.t2_array = [Aniso_19.t2]./1000;
% aRFoptions.w1_in = w1; 
% aRFoptions.w3_in = w3; %using Aniso data.w3 requires transpose
aRFoptions.useParallel=false;  
aRFoptions.nboot = 100;

% wobbling needs lsfOptions and aRFoptions 
% one angles and times are coming from a global fit of the anisotropy data
% at the center of the band
lsfParams(1) = fitParamBnd('Delta_cm',      17,     0.01,  50,     ''); %total linewidth 
% lsfParams(2) = fitParamBnd('theta0_deg',  24.7,     0,     180,    ''); %inertial cone   
lsfParams(3) = fitParamBnd('tr1',         34.4,     1,     100,    ''); %from fitting aniso
lsfParams(5) = fitParamBnd('T2'        ,  1.20,     1,     50,     ''); %from averaging T1 times
lsfParams(6) = fitParamBnd('ampSSD'    ,  0.90,     0,     1 ,     ''); %from previous section of code (one cone)
lsfParams(7) = fitParamBnd('tauSSD'    ,    10,     1,    950,     ''); %from previous section of code (one cone)

% best fit from bootstrapping CLS in clintonanalysis.m using
% RISDFitWobbling1inertial1cone1SSD.m . SSD terms come near the end as a
% CLS correction

aRFoptions.pol = 'para';
%aRFoptions.pol = 'perp';
aRFoptions.order = 1;

% the lsf needs params from the aRFoptions also to set up the numerical integration time
lsf = lsfRISDwobblingHL1diff1SSDNIBnd(lsfParams,'free',aRFoptions);

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
s(1).dataMatrix = prepareGlobalFitData(crop_data(1,:))/(2e2);

%below should be for perp
aRFoptions.pol = 'perp';
lsf = lsfRISDwobblingHL1diff1SSDNIBnd(lsfParams,'free',aRFoptions); %have to reinitialize this so that the 'perp' get applied
fun_array = {oPerp}; 
dynOther = additionalDynamics(label,fun_array,ind_array);%has no free parameters (gets them from lsf)
aRFoptions.damping = lsf;
aRFoptions.dynOther = dynOther;
s(2) = aRFCO2Bnd(aRFoptions); %make a second s for perp spectra
s(2).dataMatrix = prepareGlobalFitData(crop_data(2,:))/(2e2);

t=linspace(0,15,1000);
p = lsf.copyParamValuesToParamStruct;
figure(31),plot(t,lsf.R.para(t,p)),title('FFCF R para and perp')
hold on
plot(t,lsf.R.perp(t,p))
hold off
figure(32),
for ii = 1:4
    subplot(4,1,ii)
    plot(t,lsf.L_l{ii}(t,p)),set(gca,'Ylim',[0,1]),title(['OCF L_{',num2str(ii),'}'])
end
% figure(500)
% hold on
% plot(t,0.4.*lsf.L_l{2}(t,p),'k')
% hold off
% plot([Aniso_19.t2]./1000,[Aniso_19.area],'x')

m = multiMeasurement(s);

%% simulate via multiMeasurement

%simulated
tic
m.s(1) = m.s(1).calcSpectrum(m.s(1).p0);
% matrixPlot2DIR(s(1).simMatrix,s(1).w1_in,s(1).w3_in,s(1).t2_array,[5 6],'height_px',100,'zlimit',1) %sim

%perpendicular
m.s(2) = m.s(2).calcSpectrum(m.s(2).p0);
% matrixPlot2DIR(s(2).simMatrix,s(2).w1_in,s(2).w3_in,s(2).t2_array,[5 6],'height_px',100,'zlimit',1) %sim
toc
%% chi by eye?
matrixPlot2DIR(m.s(1).simMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit',1); %sim
matrixPlot2DIR(m.s(2).simMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit',1); %sim
% 
matrixPlot2DIR(m.s(1).dataMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit',1); %data
% matrixPlot2DIR(m.s(2).dataMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit',1); %data

%% minimize error?
% m.globalFit %takes waaaaay too long to even start up the calculations?

%% output
% % 
% matrixPlot2DIR(m.s(1).simMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit'); %sim
% matrixPlot2DIR(m.s(2).simMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit'); %sim

% matrixPlot2DIR(m.s(1).dataMatrix,m.s(1).w1_in,m.s(1).w3_in,m.s(1).t2_array,[5 6],'height_px',100,'zlimit'); %real
% matrixPlot2DIR(m.s(2).dataMatrix,m.s(2).w1_in,m.s(2).w3_in,m.s(2).t2_array,[5 6],'height_px',100,'zlimit'); %real
%cool this thing works for once You just need to z-limit the spectra to
%see that they should have the same shape



%% time to get CLS values from the simMatrices?

zzzz_data = struct('w1',num2cell(1:28), 'w3',num2cell(1:28),'R',num2cell(1:28),'t2',num2cell(1:28));

for ii = 1 : 28
zzzz_data(ii).R = m.s(1).simMatrix(:,:,ii);
zzzz_data(ii).w1 = m.s(1).w1_in;
zzzz_data(ii).w3 = m.s(1).w3_in;
zzzz_data(ii).t2 = m.s(1).t2_array(ii).*1000;
end

d_zzzz = cropData(zzzz_data(1:end),w1,w3);
t2_array_zzzz =  [d_zzzz.t2]./1000;

zzxx_data = struct('w1',num2cell(1:28), 'w3',num2cell(1:28),'R',num2cell(1:28),'t2',num2cell(1:28));

for ii = 1 : 28
zzxx_data(ii).R = m.s(2).simMatrix(:,:,ii);
zzxx_data(ii).w1 = m.s(2).w1_in;
zzxx_data(ii).w3 = m.s(2).w3_in;
zzxx_data(ii).t2 = m.s(2).t2_array(ii).*1000;
end
  %% let's try and figure out how the ESA peak intensity changes as a function of time (normalized?)
  
%   for ii = 1 : length(t2_array_zzzz)
%       I_red_spara(ii) = zzzz_data(ii).R(6,22)./zzzz_data(1).R(6,22);
%       I_red_sperp(ii) = zzxx_data(ii).R(6,22)./zzxx_data(1).R(6,22);
%   
%       I_red_epara(ii) = crop_data(1,ii).R(6,21)./crop_data(1,1).R(6,21);
%       I_red_eperp(ii) = crop_data(2,ii).R(6,21)./crop_data(2,1).R(6,21);
%   end
%   
%   figure(501)
%   plot(t2_array_zzzz,I_red_spara,'ro')
%   hold on
%   plot(t2_array_zzzz,I_red_sperp,'rx')
%   plot(t2_array_zzzz,I_red_epara,'bo')
%   plot(t2_array_zzzz,I_red_eperp,'bx')
%       xlabel('t_2 (ps)')
%       ylabel('Normalized Intensity')
%       legend({'Sim para','Sim perp','Exp para','Exp perp'})
%       title('ESA Intensity')
%       xlim([0 210])
%       ylim([0 1.2])
%   hold off
  
  %% let's try and figure out how the GSB/SE peak intensity changes as a function of time (normalized?)
  
  for ii = 1 : length(t2_array_zzzz)
      I_blue_spara_norm(ii) = zzzz_data(ii).R(30,18)./zzzz_data(1).R(30,18);
      I_blue_sperp_norm(ii) = zzxx_data(ii).R(30,18)./zzxx_data(1).R(30,18);

%       I_blue_epara_norm(ii) = crop_data(1,ii).R(20,22)./crop_data(1,1).R(20,22);
%       I_blue_eperp_norm(ii) = crop_data(2,ii).R(20,22)./crop_data(2,1).R(20,22);
  end
  
  figure(301)
  plot(t2_array_zzzz,I_blue_spara_norm,'ro')
  hold on
  plot(t2_array_zzzz,I_blue_sperp_norm,'rx')
%   plot(t2_array_zzzz,I_blue_epara_norm,'bo')
%   plot(t2_array_zzzz,I_blue_eperp_norm,'bx')
      xlabel('t_2 (ps)')
      ylabel('Normalized Intensity')
      legend({'Sim para','Sim perp'})
      title('GSB/SE Intensity')
      xlim([0 12])
      ylim([0 1.2])
  hold off

  %% check the anisotropy and isotropic decays
% p=struct();

% p.k_u     = dynParams1(1).value;
% p.k_hgs10 = dynParams1(2).value;
% p.k_r     = dynParams1(3).value;
% p.dE_cm = 590;
% p.temperature = 290 ;
% p.Delta_cm   = lsfParams(1).value; %total linewidth
% % p.theta0_deg = lsfParams(2).value; %inertial cone   
% p.tr1        = lsfParams(3).value; %cone 1 time
% p.theta1_deg = lsfParams(4).value; %cone 1 angle
% p.tr2        = lsfParams(5).value; %cone 2
% p.theta2_deg = lsfParams(6).value; %cone 2 angle
% p.tr3        = lsfParams(7).value; %diffusive time or t3
% 
%   for ii = 1 : length(t2_array_zzzz)
%       I_blue_spara_unnorm(ii) = zzzz_data(ii).R(20,22);
%       I_blue_sperp_unnorm(ii) = zzxx_data(ii).R(20,22);
%       
%       I_blue_epara_unnorm(ii) = crop_data(1,ii).R(20,22);
%       I_blue_eperp_unnorm(ii) = crop_data(2,ii).R(20,22);
%   end
   
% 
% figure(2001)
% semilogx(t2_array_zzzz,(I_blue_epara_unnorm-I_blue_eperp_unnorm)./(I_blue_epara_unnorm+2.*I_blue_eperp_unnorm),'x')
% hold on
% semilogx(t2_array_zzzz,(I_blue_spara_unnorm-I_blue_sperp_unnorm)./(I_blue_spara_unnorm+2.*I_blue_sperp_unnorm),'o')
% semilogx(t2_array_zzzz,0.4.*L_l{2}(t2_array_zzzz,p))
%     xlabel('t_2 (ps)')
%     ylabel('anisotropy')
%     legend('Exp','Sim','L_2')
% hold off


% figure(2002)
% % semilogx(t2_array_zzzz,(I_blue_epara_norm+2.*I_blue_eperp_norm)./3,'rx')
% 
% semilogx(t2_array_zzzz,(lsf.R.para(t2_array_zzzz,p)+2.*lsf.R.perp(t2_array_zzzz,p))./3,'r')
% hold on
% semilogx(t2_array_zzzz,(I_blue_epara_norm+2.*I_blue_eperp_norm)./3,'ro')
% semilogx(t2_array_zzzz,(I_blue_spara_norm+2.*I_blue_sperp_norm)./3,'bo')
% semilogx(t2_array_zzzz,L_l{1}(t2_array_zzzz,p),'k')
%     xlabel('t_2 (ps)')
%     ylabel('isotropic decay')
% %     legend('Rpara+2Rperp','L_1')
% hold off


%%

% HWHM = 21;
% pfOptions.range1 = [2059 - HWHM   2059 + HWHM];
pfOptions.range1 = [2040 2082]; %range from paper1

pfOptions.range3 = [1990 2100]; 
view2DIRdata(zzzz_data(1),pfOptions.range1,pfOptions.range3);
pfOptions.fitfcn = fittype('- a1.*voigt(w,center,wg,wl) + a2.*voigt(w,center-anh,wg,wl) +c', ...
    'coeff',{'center','a1','a2','wg','wl','anh','c'},'indep',{'w'});
pfOptions.startpoint = [2059   0.1    0.1     7    7    15   0    ];
pfOptions.lb = [        2030   0      0       0    0    10   0    ];
pfOptions.ub = [        2100   Inf    Inf     100  100  30   1000 ];
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
% HWHM = 21;
% pfOptions.range1 = [2059 - HWHM   2059 + HWHM];
pfOptions.range1 = [2040 2082]; %range from paper1

pfOptions.range3 = [1990 2100];  
view2DIRdata(zzzz_data(1),pfOptions.range1,pfOptions.range3);
pfOptions.fitfcn = fittype('- a1.*voigt(w,center,wg,wl) + a2.*voigt(w,center-anh,wg,wl) + c', ...
    'coeff',{'center','a1','a2','wg','wl','anh','c'},'indep',{'w'});
pfOptions.startpoint = [2059   0.1    0.1     5    5    15   0    ];
pfOptions.lb = [        2030   0      0       0    0    10   0    ];
pfOptions.ub = [        2100   Inf    Inf     100  100  40   1000 ];
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
% 
% cfOptions(1).startpoint = [0.3  3    0.0 ];
% cfOptions(1).lb =         [0    0.1  0    ];
% cfOptions(1).ub =         [1    500  1    ];
% cfOptions(1).fitfcn = fittype(@(a1,t1,c,x) ...
%     a1.*exp(-x./t1) + c, ...
%     'coeff', {'a1', 't1', 'c'}, 'indep', {'x'});
% cfOptions(1).flag_plot = 0; %0=off 1=on  
% 
cfOptions(2).startpoint = [0.3 0.3 3    100 0.0 ];
cfOptions(2).lb =         [0   0   0.1  10  0    ];
cfOptions(2).ub =         [1   1   10   800 1    ];
cfOptions(2).fitfcn = fittype(@(a1,a2,t1,t2,c,x) ...
    a1.*exp(-x./t1) +  a2*exp(-x./t2) + c, ...
    'coeff', {'a1', 'a2', 't1','t2', 'c'}, 'indep', {'x'});
cfOptions(2).flag_plot = 0; %0=off 1=on 

cfOptions(3).startpoint = [0.2 0.2 0.2 3    10  100 ];
cfOptions(3).lb =         [0   0   0   0.1  5   50  ];
cfOptions(3).ub =         [1   1   1   5    50  Inf];
cfOptions(3).fitfcn = fittype(@(a1,a2,a3,t1,t2,t3,x) ...
      a1.*exp(-x./t1) +  a2*exp(-x./t2)+a3.*exp(-x./t3), ...
      'coeff', {'a1', 'a2','a3', 't1','t2','t3'}, 'indep', {'x'});
cfOptions(3).flag_plot = 0; %0=off 1=on

for ii = 2:3
    corr_struct_zzzz = corrFcnFit(t2_array_zzzz,c2_zzzz,c2_std_zzzz,cfOptions(ii));
    CFFit_zzzz(ii) = corr_struct_zzzz;
end
for ii = 2:3
    corr_struct_zzxx = corrFcnFit(t2_array_zzxx,c2_zzxx,c2_std_zzxx,cfOptions(ii));
    CFFit_zzxx(ii) = corr_struct_zzxx;
end

%%

figure(302) %biexponential
  set(gcf,'color','w')
%   plot(CFFit_zzzz(3).t2(:),(CFFit_zzzz(3).fitresult(CFFit_zzzz(3).t2(:)))); 
  hold on
%   errorbar(CFFit_zzzz(2).t2(:),CFFit_zzzz(2).c2,CFFit_zzzz(2).c2_std,'bo'); 
  
%   plot(CFFit_zzxx(3).t2(:),(CFFit_zzxx(3).fitresult(CFFit_zzxx(3).t2(:)))); 
%   errorbar(CFFit_zzxx(2).t2(:),CFFit_zzxx(2).c2,CFFit_zzxx(2).c2_std,'bx');
%   plot(t2_array_zzzz,lsf.R.para(t2_array_zzzz,p))
%   plot(t2_array_zzzz,lsf.R.perp(t2_array_zzzz,p))
plot(CFFit_CLS_LIFE006_para_19C.t2(:),CFFit_CLS_LIFE006_para_19C.fitresult(CFFit_CLS_LIFE006_para_19C.t2(:)),'k-');
%   errorbar(Bmimn2oRTv2avesepcal.analysis.zzzz(3).t2(:),Bmimn2oRTv2avesepcal.analysis.zzzz(3).c2(:),Bmimn2oRTv2avesepcal.analysis.zzzz(3).c2_std(:),'bo');
%   
  plot(CFFit_CLS_LIFE006_perp_19C.t2(:),CFFit_CLS_LIFE006_perp_19C.fitresult(CFFit_CLS_LIFE006_perp_19C.t2(:)),'k--');
%   errorbar(Bmimn2oRTv2avesepcal.analysis.zzxx(3).t2(:),Bmimn2oRTv2avesepcal.analysis.zzxx(3).c2(:),Bmimn2oRTv2avesepcal.analysis.zzxx(3).c2_std(:),'ro');
%  
%SSD on data points
% uh = (0.901.*exp(-(CFFit_zzzz(2).t2)./29.5)).*CFFit_zzzz(2).c2; %values for amp and tauSSD from bootstrapping
% uhhh = (0.901.*exp(-(CFFit_zzzz(2).t2)./29.5)).*CFFit_zzxx(2).c2; %CLS values in clintonanalysis.m
rpara = lsf.R.para(t2_array_zzzz,p);
rperp = lsf.R.perp(t2_array_zzzz,p);
%SSD on R.para and Rperp
uh = (0.8.*exp(-(t2_array_zzzz)./15)).*rpara; %values for amp and tauSSD from bootstrapping
uhhh = (0.8.*exp(-(t2_array_zzzz)./15)).*rperp; %CLS values in clintonanalysis.m

plot(CFFit_zzzz(3).t2,rpara,'r-')
plot(CFFit_zzxx(3).t2,rperp,'r--')
plot(CFFit_zzzz(3).t2,uh,'ro')
plot(CFFit_zzxx(3).t2,uhhh,'rx')
xlabel('t_2 (ps)')
ylabel('c_2')
hold off
%% do some anisotropy

% for ii = 1:length(t2_array_zzzz)
% paramax(ii) = Bmimn2oRTv2avesepcal.g.para(ii).R(6,9);
% perpmax(ii) = Bmimn2oRTv2avesepcal.g.perp(ii).R(6,9); %all for real data
% 
% paracenter(ii) = Bmimn2oRTv2avesepcal.g.para(ii).R(5,7);
% perpcenter(ii) = Bmimn2oRTv2avesepcal.g.perp(ii).R(5,7);
% end
% 
% anisomax = (paramax - perpmax)./(paramax + 2.*perpmax);
% anisopeak = (paracenter - perpcenter)./(paracenter + 2.*perpcenter);
% t2 = t2_array_zzzz;

%% anisotropy of sim data
w1a = [2058 2062];
w3a = [2058 2062];
simpara = (prepareGlobalFitData(cropData(zzzz_data(1:end),w1a,w3a)));
simperp = (prepareGlobalFitData(cropData(zzxx_data(1:end),w1a,w3a)));

simaniso = (simpara - simperp)./(simpara + 2.*simperp);
simanisopeak = squeeze(simaniso(1,1,:))';
% figure(300)
%         [C,h]=contourf(cropData(zzzz_data(1),w1,w3).w1,cropData(zzzz_data(1),w1,w3).w3,simaniso(:,:,1),(0:0.02:0.7));
%         set(gcf,'PaperUnits','centimeters','PaperPosition',[5.0 5.0 7.5 5.0]);
%         title(sprintf('%0.1f ps',t2(1)));
%         xlabel('\omega_1/2\pic (cm^{-1})');
%         ylabel('\omega_3/2\pic (cm^{-1})');
%         clabel(C,h);

%%

figure(303);clf
plot([Aniso_19.t2]./1000,[Aniso_19.area],'r*')
hold on
% plot(t2_array_zzzz,simanisopeak,'bo')
plot(t2_array_zzzz,0.4.*lsf.L_l{2}(t2_array_zzzz,p),'k-')
xlim([0 12])
ylim([0 0.41])
xlabel('t_2 (ps)')
ylabel('Anisotropy')
box off
hold off
% figure(300)
%         [C,h]=contourf(cropData(zzzz_data(1),w1,w3).w1,cropData(zzzz_data(1),w1,w3).w3,simaniso(:,:,1),(0:0.02:0.7));
%         set(gcf,'PaperUnits','centimeters','PaperPosition',[5.0 5.0 7.5 5.0]);
%         title(sprintf('%0.1f ps',t2(1)));
%         xlabel('\omega_1/2\pic (cm^{-1})');
%         ylabel('\omega_3/2\pic (cm^{-1})');
%         clabel(C,h);
