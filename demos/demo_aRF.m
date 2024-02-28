%% demo of new OO-spectrum functionality
% this set of tools allows the calculation of 2D-IR spectra based on
% response functions using the truncation of the cumulant expansion at
% second order. Functionality is provided to fit the response functions to
% experimental data. There are several pieces: 
%
% fit parameters
% --------------
% 
% All parameters that can be optimized to fit the spectra to the data are
% defined through the fitParam class for open optimization or fitParamBnd
% for bounded (upper/lower limits) optimization.
%
% lineshape function g(t)
% -----------------------
%
% The lineshape function classes allow the definition of the parameters and
% functional form of the lineshape function, which determines the peak
% shape and spectral diffusion. These functions are all labeled by the
% number and kind of terms in the FFCF. The classes come in two flavors
%
% - lsf1exp: 1 exponential (Kubo) FFCF
% - lsf1exp1fast: 1 homogeneous mode (T2) and 1 exponential (Kubo) FFCF
% - lsf2exp1fast: biexponential plus homogenous (T2)
% - lsf3exp1fast: triexponential plus homogenous (T2)
%
% The parameters for the lsf are set by an lsfParams structure.
%
% analytical response function options
% ------------------------------------
%
% There are two groups or kinds of options. First, are the options that are
% constants for the calcuation, i.e. not fit parameters, e.g. things like
% the time axes and amount of zero-padding. Second, there are the fit
% parameters. The fit parameters are all specified using the fitParam
% class.
%
% analytical response function model
% ----------------------------------
%
% These models include different levels of complexity in a 2D spectrum.
% They run from the most simple, like a two-level system, to very
% complicated with many peaks in the spectrum. Each model provides
% different functionality. Most models come in two 'flavors' -- one which
% uses fminsearch to find the optimal fit and parameters are unbounded, and
% one which uses fminbnd and uses parameters with explicity upper and lower
% limits (tends to be the most useful). All initial parameters are listed
% in the field p0. After fitting, the best fit parameters are in field
% pfit. 
%
% Brief list of coded aRFs:
%
% - aRF: this abstract class defines the calling interface that all models
% must follow, but it does not calculate a spectrum.
%
% - aRFTwoLevelSystem, aRFTLSBnd: Only a single GSB peak calculated from
% rephasing and non-rephasing diagrams. This is the simplest model.
% 
% - aRFWAO, aRFWAOBnd: A model of a weakly anharmonic oscillator, there are
% now the classic six diagrams with rephasing and non-rephasing for GSB,
% SE, and ESA. There are two additional diagrams for population relaxation.
%
% - aRFCO2, aRFCO2Bnd: a model of the CO2 spectrum with the main band, the
% shoulder from manifold with 1-quantum of bend, cross-peaks due to thermal
% excitation of the bend, and direct and cross population relaxation.
% 
% - aRFEAN, aRFEANBnd: a model of the SCN- in EAN spectrum with two
% overlapping subensembles of SCN.
%
% additional dynamics
% -------------------
%
% These additional dynamics are specified to multiply each of the Feyman
% diagrams. This can treat orientational dynamics or population transfer
% during t2. 
%
% - additionalDynamics: creates functions that multiply the diagrams but
% rely on parameters defined in the model.
%
% - aDArrayBnd: creates functions that multiply the diagrams but
% also provide new (additional) parameters that are not defined in the
% model.
%
%
% multi measurement
% -----------------
%
% This is a wrapper class so that multiple experiments can be fit at the
% same time with the same set of parameters, for example, to polarization
% resolved experiments with parallel and perpendicular spectra.
%
% feynmanDiagram
% --------------
%
% This class provides functionality to calculate response functions for a
% particular diagram, rephasing or non-rephasing. This should be used only
% when defining a new aRF model. 
%

%% simple parameters
% here is an example of making a new fit parameter instance with named
% fields
p = fitParam % an empty parameter to see fields
p = fitParam('anAmplitude',3,'free') % a free parameter
pp = fitParam('anotherThing',90,'fixed') %a fixed parameter

%% bounded parameters
p = fitParamBnd('Delta1_cm',20,10,30,'');
pp = fitParamBnd('tau1',2,0.5,5,'');

%%
clear p pp

%% lineshapes (unbounded)
% here is an example of a simple lineshape function

lsf = lsf1exp %if declared empty returns a list of its parameters in the name field

% show the name field by itself
lsf.name

% print names
fprintf(1,'The names of the parameters for this model are: %s\n',lsf.name{:})

%% set a simple lsf unbounded
lsfParams = fitParam('Delta1_cm',20,'');
lsfParams(2) = fitParam('tau1',2,'');
lsf = lsf1exp(lsfParams,'free')

lsf.name %returns the names of the fields
lsf.value %returns the values of the parameters

%% set a simple lsf bounded
lsfParams = fitParamBnd('Delta1_cm',20,10,30,'');
lsfParams(2) = fitParamBnd('tau1',2,0.5,5,'');
lsf = lsf1expBnd(lsfParams,'free')

%% simplest example of calculating a spectrum (single t2)
clear aRFoptions lsf s

% set lsf
lsfParams = fitParam('Delta1_cm',20,'');
lsfParams(2) = fitParam('tau1',2,'');
lsf = lsf1exp(lsfParams,'free')

% set fixed params
aRFoptions.dt = 0.1; % time step
aRFoptions.n_t = 64; % number of steps
aRFoptions.n_zp = 2*aRFoptions.n_t; % zero-padding
aRFoptions.flag_rotating_frame = true; % rotating frame (should be true always)
aRFoptions.t2_array = 0; % list of t2 times to calculate
aRFoptions.w1_in = 2200:10:2400; % freq axes of interest (should match experiment)
aRFoptions.w3_in = 2200:10:2400;
%
aRFoptions.w_01_cm = fitParam('w_01_cm',2350,'fixed');
aRFoptions.mu01sq = fitParam('mu01sq',1,'fixed');
aRFoptions.phase_deg = fitParam('phase_deg',0,'free');
aRFoptions.damping = lsf;

% this initializes everything but doesn't calculate anything (yet)
s = aRFTwoLevelSystem(aRFoptions)

% the list of all parameters in the model 
disp('all parameters in the model:')
s.allFitParamNames

% the list of all parameters in the model 
disp('free parameters in the model:')
s.freeFitParamNames

%% the components of the calculation can be interrogated if needed
% one can compare these plots to the figures in H&Z

s.t2 = 0; % the t2 for the spectrum that will be calculated
s = s.makeResponseFunctions(s.paramStruct); % setup the response functions
s = s.calcDiagramsTime; %calculate each diagram in the time domain

figure(1),contourf(s.diagrams(1).R) %re
figure(2),contourf(s.diagrams(2).R) %nonre

s = s.calcDiagramsFreq(s.n_zp); % calculate the frequency domain

figure(3),clf,contourf(s.diagrams(1).R) %re
figure(4),clf,contourf(s.diagrams(2).R) %nonre

%% the total spectrum can be calculated for the initial parameters
s = s.calcSpectrum(s.p0)
figure(5),clf,my2dPlot(s.W1,s.W3,s.spec)

%% the total spectrum can be calculated in one go for an arbitraty set of parameters p 
p = [0 25 2]; % have to match the order listed in s.freeFitParamNames if you do it this way
s = s.calcSpectrum(p)
figure(6),clf,my2dPlot(s.W1,s.W3,s.spec)

%% let's try an array of t2 values
% now the individual spectra are cropped to the roi and placed in the
% matrix s.simMatrix

% set lsf
lsfParams = fitParam('Delta1_cm',20,'');
lsfParams(2) = fitParam('tau1',2,'');
lsf = lsf1exp(lsfParams,'free')

% fixed
aRFoptions.dt = 0.1;
aRFoptions.n_t = 64;
aRFoptions.n_zp = 2*aRFoptions.n_t;
aRFoptions.flag_rotating_frame = true;
aRFoptions.t2_array = [0 10 100];
aRFoptions.w1_in = 2250:10:2450;
aRFoptions.w3_in = 2250:10:2450;
% params
aRFoptions.w_01_cm = fitParam('w_01_cm',2350,'fixed');
aRFoptions.mu01sq = fitParam('mu01sq',1,'fixed');
aRFoptions.phase_deg = fitParam('phase_deg',0,'free');
aRFoptions.damping = lsf;

% initialize
s = aRFTwoLevelSystem(aRFoptions);

% calculate
s = s.calcSpectrum(s.p0)

% plot
i_spectrum = 1; % pick which spectrum to look at
figure(5),clf,my2dPlot(s.w1_in,s.w3_in,squeeze(s.simMatrix(:,:,i_spectrum))) 

%% the same structure works with bounded parameters (TLS)
clear aRFoptions lsf lsfParams s

% lsf
lsfParams = fitParamBnd('Delta1_cm',20,10,30,'');
lsfParams(2) = fitParamBnd('tau1',2,0.5,5,'');
lsf = lsf1expBnd(lsfParams,'free')

%fixed
aRFoptions.dt = 0.1;
aRFoptions.n_t = 64;
aRFoptions.n_zp = 2*aRFoptions.n_t;
aRFoptions.flag_rotating_frame = true;
aRFoptions.t2_array = [0 10 100];
aRFoptions.w1_in = 2250:10:2450;
aRFoptions.w3_in = 2250:10:2450;
% fit params
aRFoptions.w_01_cm = fitParamBnd('w_01_cm',2350,2310,2400,'fixed');
aRFoptions.mu01sq = fitParamBnd('mu01sq',2,0,inf,'free');
aRFoptions.phase_deg = fitParamBnd('phase_deg',0,-90,90,'free');
aRFoptions.damping = lsf;

% initialize
s = aRFTLSBnd(aRFoptions)

% calc initial spectrum
s = s.calcSpectrum(s.p0)

% plot
figure(5),clf,my2dPlot(s.w1_in,s.w3_in,squeeze(s.simMatrix(:,:,1)),'pumpprobe',false) 
 
%% bounded WAO
clear aRFoptions lsf lsfParams s

% lsf
lsfParams = fitParamBnd('Delta1_cm',20,10,30,'');
lsfParams(2) = fitParamBnd('tau1',2,0.5,5,'');
lsf = lsf1expBnd(lsfParams,'free')

% fixed
aRFoptions.dt = 0.1;
aRFoptions.n_t = 64;
aRFoptions.n_zp = 2*aRFoptions.n_t;
aRFoptions.flag_rotating_frame = true;
aRFoptions.t2_array = [0 10 100];
aRFoptions.w1_in = 2250:10:2450;
aRFoptions.w3_in = 2250:10:2450;
% fit params
aRFoptions.w_01_cm = fitParamBnd('w_01_cm',2350,2310,2400,'fixed');
aRFoptions.mu01sq = fitParamBnd('mu01sq',1000,0,inf,'free');
aRFoptions.mu12sqRatio = fitParamBnd('mu12sqRatio',1,0,inf,'free');
aRFoptions.anh_cm = fitParamBnd('anh_cm',50,0,inf,'free');
aRFoptions.phase_deg = fitParamBnd('phase_deg',0,-90,90,'free');
aRFoptions.damping = lsf;

% initialize
s = aRFWAOBnd(aRFoptions);

% calculate initial spectrum
s = s.calcSpectrum(s.p0)

% plot
figure(5),clf,my2dPlot(s.w1_in,s.w3_in,squeeze(s.simMatrix(:,:,1)),'pumpprobe',false) 
 
%% um
s = aRFTLSBnd(aRFoptions)
s = s.calcSpectrum([4 -5 12 1])
figure(4),clf,my2dPlot(s.w1_in,s.w3_in,squeeze(s.simMatrix(:,:,1)),'pumpprobe',false) 
s.dataMatrix = s.simMatrix;
%s.err_fun(s.p0);
s = s.globalFit;
s = s.calcSpectrum(s.pfit);
figure(5),clf,my2dPlot(s.w1_in,s.w3_in,squeeze(s.simMatrix(:,:,1)),'pumpprobe',false) 


%% start working on fitting
% now we can look at how to fit spectra to mock experimental data just to
% demonstrate

%% bounded WAO
clear aRFoptions lsf lsfParams s

% lsf
lsfParams = fitParamBnd('Delta1_cm',20,10,30,'');
lsfParams(2) = fitParamBnd('tau1',2,0.5,5,'');
lsf = lsf1expBnd(lsfParams,'free')

% fixed
aRFoptions.dt = 0.1;
aRFoptions.n_t = 64;
aRFoptions.n_zp = 2*aRFoptions.n_t;
aRFoptions.flag_rotating_frame = true;
aRFoptions.t2_array = [0 10 100];
aRFoptions.w1_in = 2250:10:2450;
aRFoptions.w3_in = 2250:10:2450;
% fit params
aRFoptions.w_01_cm = fitParamBnd('w_01_cm',2350,2310,2400,'fixed');
aRFoptions.mu01sq = fitParamBnd('mu01sq',1000,0,inf,'free');
aRFoptions.mu12sqRatio = fitParamBnd('mu12sqRatio',1,0,inf,'free');
aRFoptions.anh_cm = fitParamBnd('anh_cm',50,0,inf,'free');
aRFoptions.phase_deg = fitParamBnd('phase_deg',0,-90,90,'free');
aRFoptions.damping = lsf;

% initialize
s = aRFWAOBnd(aRFoptions);

% calculate initial spectrum
s = s.calcSpectrum(s.p0)

% plot
figure(5),clf,my2dPlot(s.w1_in,s.w3_in,squeeze(s.simMatrix(:,:,1)),'pumpprobe',false) 
 
% calculate a spectrum from some different parameters than what we start
% with (not p0 but close)
ptarget = [900 1.5 35 0 12 1];
s = s.calcSpectrum(ptarget);

% plot the target spectrum
figure(4),clf,my2dPlot(s.w1_in,s.w3_in,squeeze(s.simMatrix(:,:,1)),'pumpprobe',false) 
drawnow

% copy the simulated spectrum to the "dataMatrix" so that this is what the
% global fit will adjust the parameters to reproduce
s.dataMatrix = s.simMatrix;

% do the fit
s = s.globalFit;

% look at the resulting fit
s = s.calcSpectrum(s.pfit);
figure(5),clf,my2dPlot(s.w1_in,s.w3_in,squeeze(s.simMatrix(:,:,1)),'pumpprobe',false) 

% show target values
disp('Target values')
s.prettyPrintFreeParamValues(ptarget)

% show the fit result values
disp('Fit values')
s.prettyPrintFreeParamValues(s.pfit)

%% now we can add some complexity by adding population and orientational
% relaxation processes with the additionalDynamics class
% use the aDArrayBnd class that has parameters packaged

%good idea to clear this so fields don't cascade
clear aRFoptions lsf lsfParams dyn dynParams dynParams2 s 

% new additional dynamics
dyn = aDArrayBnd;%first time call empty, later can call with indices

%T1 population relaxation
label = 'T1_rel';
f1 = @(T1,t2,T3,p) exp(-T1./(2*p.T1_rel2)).*exp(-T3./(2*p.T1_rel2));
f2 = @(T1,t2,T3,p) exp(-T1./(2*p.T1_rel2)).* exp(-t2./p.T1_rel2).*exp(-T3./(2*p.T1_rel2));
f3 = @(T1,t2,T3,p) exp(-T1./(2*p.T1_rel2)).*(1-exp(-t2./p.T1_rel2)).*exp(-T3./(2*p.T1_rel2));
fun_array = {f1, f2, f3};
ind_array = {[1 4], [2 3 5 6], [7 8]};
dynParams1    = fitParamBnd('T1_rel2',10,1,100,'');
dyn(1) = aDArrayBnd('T1_rel',fun_array,ind_array,dynParams1,'free');

%Orientation relaxation
label = 'Orientation_rel';
oPara = @(T1,t2,T3,p) (1/9).*exp(-T1./(3.*p.tau_o2)).*(1+4/5.*exp(-t2./p.tau_o2)).*exp(-T3./(3*p.tau_o2));
oPerp = @(T1,t2,T3,p) (1/9).*exp(-T1./(3.*p.tau_o2)).*(1-2/5.*exp(-t2./p.tau_o2)).*exp(-T3./(3*p.tau_o2));
oCrossed = @(T1,t2,T3,p) 1/15.*exp(-T1./(3.*p.tau_o2)).*exp(-t2./p.tau_o2).*exp(-T3./(3*p.tau_o2));
oMagic = @(T1,t2,T3,p) 1/9.*exp(-T1./(3.*p.tau_o2)).*exp(-T3./(3*p.tau_o2));%need to checkprefactors
fun_array = {oPara};
ind_array = {[1:8]};
dynParams2 = fitParamBnd('tau_o2',10,1,100,'');
dyn(2) = aDArrayBnd(label,fun_array,ind_array,dynParams2,'fixed');


lsfParams = fitParamBnd('Delta1_cm',20,10,30,'');
lsfParams(2) = fitParamBnd('tau1',2,0.5,5,'');
lsf = lsf1expBnd(lsfParams,'free');

clear aRFoptions ;
aRFoptions.dt = 0.1;
aRFoptions.n_t = 64;
aRFoptions.n_zp = 2*aRFoptions.n_t;
aRFoptions.flag_rotating_frame = true;
aRFoptions.t2_array = [0 10 100];
aRFoptions.w1_in = 2250:10:2450;
aRFoptions.w3_in = 2250:10:2450;
aRFoptions.useParallel = true;
aRFoptions.w_01_cm = fitParamBnd('w_01_cm',2350,2310,2400,'fixed');
aRFoptions.mu01sq = fitParamBnd('mu01sq',1000,0,inf,'free');
aRFoptions.mu12sqRatio = fitParamBnd('mu12sqRatio',2,0,inf,'free');
aRFoptions.anh_cm = fitParamBnd('anh_cm',50,0,inf,'free');
aRFoptions.phase_deg = fitParamBnd('phase_deg',0,-90,90,'free');
aRFoptions.damping = lsf;
aRFoptions.dyn = dyn;

% initialize
s = aRFWAOBnd(aRFoptions)

% calc spectrum
s = s.calcSpectrum(s.p0)

% plot
figure(5),clf,my2dPlot(s.w1_in,s.w3_in,squeeze(s.simMatrix(:,:,1)),'pumpprobe',false) 

%% try fitting to test
s = s.calcSpectrum(s.p0.*1.05);
figure(4),clf,my2dPlot(s.w1_in,s.w3_in,squeeze(s.simMatrix(:,:,1)),'pumpprobe',false) 
drawnow
s.dataMatrix = s.simMatrix;
%s.err_fun(s.p0);
s = s.globalFit;
s = s.calcSpectrum(s.pfit);
figure(5),clf,my2dPlot(s.w1_in,s.w3_in,squeeze(s.simMatrix(:,:,1)),'pumpprobe',false) 

%% an example of bootstrapping
s.nboot=10  
s = s.globalFitBootstrap;



%% STOP HERE the rest is very old and may break...
%% demo multiMeasurement (this is for polarization resolved experiments)

% first set up the parallel measurement. Most of these properties will be
% copied into the perp measurement, too. Only the orientational dynamics
% will be different. All parameters are shared.
dyn = additionalDynamics;%first time call empty, later can call with indices

%T1 population relaxation
label = 'T1_rel';
f1 = @(T1,t2,T3,p) exp(-T1./(2*p.T1_rel)).*exp(-T3./(2*p.T1_rel));
f2 = @(T1,t2,T3,p) exp(-T1./(2*p.T1_rel)).* exp(-t2./p.T1_rel).*exp(-T3./(2*p.T1_rel));
f3 = @(T1,t2,T3,p) exp(-T1./(2*p.T1_rel)).*(1-exp(-t2./p.T1_rel)).*exp(-T3./(2*p.T1_rel));
fun_array = {f1, f2, f3};
ind_array = {[1 4], [2 3 5 6], [7 8]};
dyn(1) = additionalDynamics(label,fun_array,ind_array);

%Orientation relaxation
label = 'Orientation_rel';
oPara = @(T1,t2,T3,p) (1/9).*exp(-T1./(3.*p.tau_o)).*(1+4/5.*exp(-t2./p.tau_o)).*exp(-T3./(3*p.tau_o));
oPerp = @(T1,t2,T3,p) (1/9).*exp(-T1./(3.*p.tau_o)).*(1-2/5.*exp(-t2./p.tau_o)).*exp(-T3./(3*p.tau_o));
oCrossed = @(T1,t2,T3,p) 1/15.*exp(-T1./(3.*p.tau_o)).*exp(-t2./p.tau_o).*exp(-T3./(3*p.tau_o));
oMagic = @(T1,t2,T3,p) 1/9.*exp(-T1./(3.*p.tau_o)).*exp(-T3./(3*p.tau_o));%need to checkprefactors
fun_array = {oPara}; %this is the important line that sets it as parallel!!!
ind_array = {[1:8]};
dyn(2) = additionalDynamics(label,fun_array,ind_array);


lsfParams = fitParamBnd('Delta1_cm',20,10,30,'');
lsfParams(2) = fitParamBnd('tau1',2,0.5,5,'');
lsf = lsf1expBnd(lsfParams,'free');

aRFoptions.dt = 0.1;
aRFoptions.n_t = 64;
aRFoptions.n_zp = 2*aRFoptions.n_t;
aRFoptions.flag_rotating_frame = true;
aRFoptions.w_01_cm = fitParamBnd('w_01_cm',2350,2310,2400,'fixed');
aRFoptions.mu01sq = fitParamBnd('mu01sq',1000,0,inf,'free');
aRFoptions.mu12sqRatio = fitParamBnd('mu12sqRatio',2,0,inf,'free');
aRFoptions.anh_cm = fitParamBnd('anh_cm',50,0,inf,'free');
aRFoptions.phase_deg = fitParamBnd('phase_deg',0,-90,90,'free');
aRFoptions.damping = lsf;
aRFoptions.dyn = dyn;
aRFoptions.T1_rel = fitParamBnd('T1_rel',10,1,100,'free');
aRFoptions.tau_o = fitParamBnd('tau_o',11,1,100,'free');
aRFoptions.t2_array = [0 10 100];
aRFoptions.w1_in = 2250:10:2450;
aRFoptions.w3_in = 2250:10:2450;

% make the first spectrum object instance
s = aRFWAOBnd(aRFoptions);

% now make the second spectrum object
fun_array = {oPerp}; %this makes it perpendicular
dyn(2) = additionalDynamics(label,fun_array,ind_array);
aRFoptions.dyn = dyn; % save as the only difference

% make second spectrum object instance
s(2) = aRFWAOBnd(aRFoptions); 

m = multiMeasurement(s);

%% try fitting
%first set up the "data" matrices (this will be assigning real data
%eventually)
p_target = [900 1.8 45 12 15 2 22 2.2];
%parallel
m.s(1) = m.s(1).calcSpectrum(p_target);
m.s(1).dataMatrix = m.s(1).simMatrix;
%perpendicular
m.s(2) = m.s(2).calcSpectrum(p_target);
m.s(2).dataMatrix = m.s(2).simMatrix;

m.globalFit




%% start CO2
dyn = additionalDynamics;%first time call empty, later can call with indices


%Orientation relaxation
label = 'Orientation_rel';
oPara = @(T1,t2,T3,p) (1/9).*exp(-T1./(3.*p.tau_o)).*(1+4/5.*exp(-t2./p.tau_o)).*exp(-T3./(3*p.tau_o));
oPerp = @(T1,t2,T3,p) (1/9).*exp(-T1./(3.*p.tau_o)).*(1-2/5.*exp(-t2./p.tau_o)).*exp(-T3./(3*p.tau_o));
oCrossed = @(T1,t2,T3,p) 1/15.*exp(-T1./(3.*p.tau_o)).*exp(-t2./p.tau_o).*exp(-T3./(3*p.tau_o));
oMagic = @(T1,t2,T3,p) 1/9.*exp(-T1./(3.*p.tau_o)).*exp(-T3./(3*p.tau_o));%need to checkprefactors
fun_array = {oPara};
ind_array = {[1:32]};
dyn(1) = additionalDynamics(label,fun_array,ind_array);

%cheat to turn some diagrams off and on
label = 'cheat';
on = @(T1,t2,T3,p) 1;
off = @(T1,t2,T3,p) 0;
fun_array = {off,on};
ind_array = {[1:6 7:12 13:18 19:24 ],[25:26  27:28 29:30 31:32]};
dyn(2) = additionalDynamics(label,fun_array,ind_array);

lsfParams = fitParamBnd('Delta1_cm',3.77,0.1,5,''); %%%check how calcd!!!
lsfParams(2) = fitParamBnd('tau1',101,10,250,'');
lsfParams(3) = fitParamBnd('T2',7.7,0.1,15,'');
lsf = lsf1exp1fastBnd(lsfParams,'free');

%
clear aRFoptions s;
aRFoptions.dt = 0.2;
aRFoptions.n_t = 64;
aRFoptions.n_zp = 2*aRFoptions.n_t;
aRFoptions.flag_rotating_frame = true;
aRFoptions.w_01_cm = fitParamBnd('w_01_cm',2337.5,2310,2400,'fixed');
aRFoptions.dw_sb_cm = fitParamBnd('dw_sb_cm',12,10,15,'fixed');
aRFoptions.mu01sq = fitParamBnd('mu01sq',1000,0,inf,'fixed');
aRFoptions.mu12sqRatio = fitParamBnd('mu12sqRatio',1.51,0,inf,'fixed');
aRFoptions.anh_cm = fitParamBnd('anh_cm',23.5,0,inf,'fixed');
aRFoptions.phase_deg = fitParamBnd('phase_deg',2,-90,90,'fixed');
aRFoptions.damping = lsf;
aRFoptions.dyn = dyn;
aRFoptions.T1_rel = fitParamBnd('T1_rel',10,1,100,'fixed');%get rid of this???
aRFoptions.tau_o = fitParamBnd('tau_o',100,1,100,'fixed');
aRFoptions.dE_cm = fitParamBnd('dE_cm',667,650,700,'fixed');
aRFoptions.temperature = fitParamBnd('temperature',296,290,305,'fixed');
aRFoptions.t2_array = [0 10 100];
aRFoptions.w1_in = 2320:0.5:2350;
aRFoptions.w3_in = 2296:2:2350;

s = aRFCO2Bnd(aRFoptions)
s = s.calcSpectrum(s.p0)
%figure(5),clf,my2dPlot(s.w1_in,s.w3_in,squeeze(s.simMatrix(:,:,1)),'pumpprobe',false) 
figure(6),clf,my2dPlot(s.w1_in,s.w3_in,squeeze(s.simMatrix(:,:,1)),'pumpprobe',false,'n_contours',20) 
matrixPlot2DIR(s.simMatrix,s.w1_in,s.w3_in,s.t2_array,[1 3])
%%
s = s.calcSpectrum([900]);
figure(4),clf,my2dPlot(s.w1_in,s.w3_in,squeeze(s.simMatrix(:,:,1)),'pumpprobe',false) 
drawnow
s.dataMatrix = s.simMatrix;
%s.err_fun(s.p0);
s = s.globalFit;
s = s.calcSpectrum(s.pfit);
figure(5),clf,my2dPlot(s.w1_in,s.w3_in,squeeze(s.simMatrix(:,:,1)),'pumpprobe',false) 

%% OLD BAD BAD
%% SCN in EAN

% first set up the parallel measurement. Most of these properties will be
% copied into the perp measurement, too. Only the orientational dynamics
% will be different. All parameters are shared.
dyn = additionalDynamics;%first time call empty, later can call with indices

%T1 population relaxation
label = 'T1_rel';
f1 = @(T1,t2,T3,p) exp(-T1./(2*p.T1_rel)).*exp(-T3./(2*p.T1_rel));
f2 = @(T1,t2,T3,p) exp(-T1./(2*p.T1_fast)).*exp(-T3./(2*p.T1_fast));
fun_array = {f1, f2};
ind_array = {[1:6], [7:12]};
dyn(1) = additionalDynamics(label,fun_array,ind_array);

%Orientation relaxation
label = 'Orientation_rel';
oPara = @(T1,t2,T3,p) (1/9).*exp(-T1./(3.*p.tau_o)).*(1+4/5.*exp(-t2./p.tau_o)).*exp(-T3./(3*p.tau_o));
oPerp = @(T1,t2,T3,p) (1/9).*exp(-T1./(3.*p.tau_o)).*(1-2/5.*exp(-t2./p.tau_o)).*exp(-T3./(3*p.tau_o));
oCrossed = @(T1,t2,T3,p) 1/15.*exp(-T1./(3.*p.tau_o)).*exp(-t2./p.tau_o).*exp(-T3./(3*p.tau_o));
oMagic = @(T1,t2,T3,p) 1/9.*exp(-T1./(3.*p.tau_o)).*exp(-T3./(3*p.tau_o));%need to checkprefactors
oPara2 = @(T1,t2,T3,p) (1/9).*exp(-T1./(3.*p.tau_o_fast)).*(1+4/5.*exp(-t2./p.tau_o_fast)).*exp(-T3./(3*p.tau_o_fast));
oPerp2 = @(T1,t2,T3,p) (1/9).*exp(-T1./(3.*p.tau_o_fast)).*(1-2/5.*exp(-t2./p.tau_o_fast)).*exp(-T3./(3*p.tau_o_fast));
oCrossed2 = @(T1,t2,T3,p) 1/15.*exp(-T1./(3.*p.tau_o_fast)).*exp(-t2./p.tau_o_fast).*exp(-T3./(3*p.tau_o_fast));
oMagic2 = @(T1,t2,T3,p) 1/9.*exp(-T1./(3.*p.tau_o_fast)).*exp(-T3./(3*p.tau_o_fast));%need to checkprefactors
fun_array = {oPara,oPara2}; %this is the important line that sets it as parallel!!!
ind_array = {[1:6],[7:12]};
dyn(2) = additionalDynamics(label,fun_array,ind_array);
%%

label = 'ffcf dynamics'
lsfParams = fitParamBnd('Delta1_cm',20,10,30,'');
lsfParams(2) = fitParamBnd('tau1',2,0.5,5,'');
lsf1 = lsf1expBnd(lsfParams,'free');
lsfParams = fitParamBnd('Delta1_cm',10,10,30,'');
lsfParams(2) = fitParamBnd('tau1',5,0.5,5,'');
lsfParams(3) = fitParamBnd('T2',1,0.1,5,'');
lsf2 = lsf1exp1fastBnd(lsfParams,'free');

%%
lsf_array = {lsf1,lsf2};
ind_array = {1:6,7:12};
%lsf = additionalDynamics(label,lsf_array,ind_array);
lsf = lsfArrayBnd(lsf_array,ind_array);


%%
 
aRFoptions.dt = 0.1;
aRFoptions.n_t = 64;
aRFoptions.n_zp = 2*aRFoptions.n_t;
aRFoptions.flag_rotating_frame = true;
aRFoptions.w_01_cm = fitParamBnd('w_01_cm',2350,2310,2400,'fixed');
aRFoptions.dw_cm = fitParamBnd('dw_cm',50,10,50,'fixed');%shift to low freq peak
aRFoptions.amp2 = fitParamBnd('amp2',1,0,1,'fixed');
aRFoptions.mu01sq = fitParamBnd('mu01sq',1000,0,inf,'free');
aRFoptions.mu12sqRatio = fitParamBnd('mu12sqRatio',2,0,inf,'free');
aRFoptions.anh_cm = fitParamBnd('anh_cm',50,0,inf,'free');
aRFoptions.phase_deg = fitParamBnd('phase_deg',0,-90,90,'free');
aRFoptions.damping = lsf;
aRFoptions.dyn = dyn;
aRFoptions.T1_rel = fitParamBnd('T1_rel',10,1,100,'free');
aRFoptions.T1_fast = fitParamBnd('T1_fast',1,1,100,'free');
aRFoptions.tau_o = fitParamBnd('tau_o',11,1,100,'free');
aRFoptions.tau_o_fast = fitParamBnd('tau_o_fast',5,1,100,'free');
aRFoptions.t2_array = [0 3 10];
aRFoptions.w1_in = 2250:10:2450;
aRFoptions.w3_in = 2250:10:2450;

% make the first spectrum object instance
s = aRFEANBnd(aRFoptions);
s = s.calcSpectrum(s.p0);
matrixPlot2DIR(s(1).simMatrix,s(1).w1_in,s(1).w3_in,s(1).t2_array,[1 3])

%%
% now make the second spectrum object
fun_array = {oPerp,oPerp2}; %this makes it perpendicular
dyn(2) = additionalDynamics(label,fun_array,ind_array);
aRFoptions.dyn = dyn; % save as the only difference

% make second spectrum object instance
s(2) = aRFEANBnd(aRFoptions); 

m = multiMeasurement(s);

% plotting
m.s(1) = m.s(1).calcSpectrum(p_target);
figure(5),clf,my2dPlot(m.s(1).w1_in,m.s(1).w3_in,squeeze(m.s(1).simMatrix(:,:,2)),'pumpprobe',false) 
drawnow
%matrixPlot2DIR(s(1).simMatrix,s(1).w1_in,s(1).w3_in,s(1).t2_array,[1 3])

%% try fitting
%first set up the "data" matrices (this will be assigning real data
%eventually)
p_target = [3 2 900 1.8 20 12 15 0 10 4];
%parallel
m.s(1) = m.s(1).calcSpectrum(p_target);
m.s(1).dataMatrix = m.s(1).simMatrix;

%perpendicular
m.s(2) = m.s(2).calcSpectrum(p_target);
m.s(2).dataMatrix = m.s(2).simMatrix;

m.globalFit

%% Try EAN again

% first set up the parallel measurement. Most of these properties will be
% copied into the perp measurement, too. Only the orientational dynamics
% will be different. All parameters are shared.
dyn = additionalDynamics;%first time call empty, later can call with indices

%T1 population relaxation
label = 'T1_rel';
f1 = @(T1,t2,T3,p) exp(-T1./(2*p.T1_1)).*exp(-t2./(p.T1_1)).*exp(-T3./(2*p.T1_1));
f2 = @(T1,t2,T3,p) exp(-T1./(2*p.T1_2)).*exp(-t2./(p.T1_2)).*exp(-T3./(2*p.T1_2));
fun_array = {f1, f2};
ind_array = {[1:6], [7:12]};
dyn(1) = additionalDynamics(label,fun_array,ind_array);

%Orientation relaxation
label = 'Orientation_rel';
oPara = @(T1,t2,T3,p) (1/9).*exp(-T1./(3.*p.tau_o_1)).*(1+4/5.*exp(-t2./p.tau_o_1)).*exp(-T3./(3*p.tau_o_1));
oPerp = @(T1,t2,T3,p) (1/9).*exp(-T1./(3.*p.tau_o_1)).*(1-2/5.*exp(-t2./p.tau_o_1)).*exp(-T3./(3*p.tau_o_1));
oCrossed = @(T1,t2,T3,p) 1/15.*exp(-T1./(3.*p.tau_o_1)).*exp(-t2./p.tau_o_1).*exp(-T3./(3*p.tau_o_1));
oMagic = @(T1,t2,T3,p) 1/9.*exp(-T1./(3.*p.tau_o_1)).*exp(-T3./(3*p.tau_o_1));%need to checkprefactors
oPara2 = @(T1,t2,T3,p) (1/9).*exp(-T1./(3.*p.tau_o_2)).*(1+4/5.*exp(-t2./p.tau_o_2)).*exp(-T3./(3*p.tau_o_2));
oPerp2 = @(T1,t2,T3,p) (1/9).*exp(-T1./(3.*p.tau_o_2)).*(1-2/5.*exp(-t2./p.tau_o_2)).*exp(-T3./(3*p.tau_o_2));
oCrossed2 = @(T1,t2,T3,p) 1/15.*exp(-T1./(3.*p.tau_o_2)).*exp(-t2./p.tau_o_2).*exp(-T3./(3*p.tau_o_2));
oMagic2 = @(T1,t2,T3,p) 1/9.*exp(-T1./(3.*p.tau_o_2)).*exp(-T3./(3*p.tau_o_2));%need to checkprefactors
fun_array = {oPara,oPara2}; %this is the important line that sets it as parallel!!!
ind_array = {[1:6],[7:12]};
dyn(2) = additionalDynamics(label,fun_array,ind_array);
%

label = 'ffcf dynamics'
lsfParams = fitParamBnd('Delta1_cm',20,10,30,'');
lsfParams(2) = fitParamBnd('tau1',15,0.5,25,'');
lsfParams(3) = fitParamBnd('T2',0.6,0.1,5,'');
lsf1 = lsf1exp1fastBnd(lsfParams,'free');
lsfParams = fitParamBnd('Delta1_cm',15,10,30,'');
lsfParams(2) = fitParamBnd('tau1',15,0.5,25,'');
lsfParams(3) = fitParamBnd('T2',0.6,0.1,5,'');
lsf2 = lsf1exp1fastBnd(lsfParams,'free');

%
lsf_array = {lsf1,lsf2};
ind_array = {1:6,7:12};
%lsf = additionalDynamics(label,lsf_array,ind_array);
lsf = lsfArrayBnd(lsf_array,ind_array,'fixed');


%
aRFoptions.dt = 0.05;
aRFoptions.n_t = 64;
aRFoptions.n_zp = 2*aRFoptions.n_t;
aRFoptions.flag_rotating_frame = true;
aRFoptions.w_01_cm = fitParamBnd('w_01_cm',2075,2060,2080,'free');
aRFoptions.dw_cm = fitParamBnd('dw_cm',50,10,50,'free');%shift to low freq peak
aRFoptions.amp2 = fitParamBnd('amp2',1,0,10,'free');
aRFoptions.mu01sq = fitParamBnd('mu01sq',1000,0,inf,'free');
aRFoptions.mu12sqRatio = fitParamBnd('mu12sqRatio',2,0,inf,'fixed');
aRFoptions.anh_cm = fitParamBnd('anh_cm',26,0,inf,'fixed');
aRFoptions.phase_deg = fitParamBnd('phase_deg',0,-90,90,'fixed');
aRFoptions.damping = lsf;
aRFoptions.dyn = dyn;
% aRFoptions.T1_1 = fitParamBnd('T1_1',5.5,1,100,'free');
% aRFoptions.T1_2 = fitParamBnd('T1_2',0.6,1,100,'free');
aRFoptions.T1_1 = fitParamBnd('T1_1',5.5,1,100,'fixed');
aRFoptions.T1_2 = fitParamBnd('T1_2',1,1,100,'fixed');
aRFoptions.tau_o_1 = fitParamBnd('tau_o_1',11,1,100,'fixed');
aRFoptions.tau_o_2 = fitParamBnd('tau_o_2',5,1,100,'fixed');
aRFoptions.t2_array = [0 1 3 10 ];
aRFoptions.w1_in = 1950:10:2150;
aRFoptions.w3_in = 1950:10:2150;

% make the first spectrum object instance
s = aRFEANBnd(aRFoptions);

%
s = s.calcSpectrum(s.p0);
matrixPlot2DIR(s(1).simMatrix,s(1).w1_in,s(1).w3_in,s(1).t2_array,[1 4])

% now make the second spectrum object
fun_array = {oPerp,oPerp2}; %this makes it perpendicular
dyn(2) = additionalDynamics(label,fun_array,ind_array);
aRFoptions.dyn = dyn; % save as the only difference

% make second spectrum object instance
s(2) = aRFEANBnd(aRFoptions); 

m = multiMeasurement(s);
%% try fitting
%first set up the "data" matrices (this will be assigning real data
%eventually)
p_target = [1.2 23 950 2073];
%parallel
m.s(1) = m.s(1).calcSpectrum(p_target);
m.s(1).dataMatrix = m.s(1).simMatrix;

%perpendicular
m.s(2) = m.s(2).calcSpectrum(p_target);
m.s(2).dataMatrix = m.s(2).simMatrix;

m = m.globalFit



