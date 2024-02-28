%% test reconstruction of spectrum from simulation
%
clear aRFoptions s lsf lsfParams dyn dynOther
lsfParams = fitParamBnd('Delta1_cm',20,10,30,'');
lsfParams(2) = fitParamBnd('tau1',2,0.5,5,'');
lsf = lsf1expBnd(lsfParams,'free');


aRFoptions.dt = 0.01;
aRFoptions.n_t = 256;
aRFoptions.n_zp = 2*aRFoptions.n_t;
aRFoptions.flag_rotating_frame = true;
aRFoptions.w_01_cm = fitParam('w_01_cm',2350,'fixed');
aRFoptions.mu01sq = fitParam('mu01sq',1,'fixed');
aRFoptions.phase_deg = fitParam('phase_deg',0,'free');
aRFoptions.damping = lsf;
aRFoptions.t2_array = 0;
aRFoptions.w1_in = 2200:10:2400;
aRFoptions.w3_in = 2200:10:2400;


s = aRFTwoLevelSystem(aRFoptions);
s = s.calcSpectrum(s.p0);

figure,clf,my2dPlot(s.w,s.w,s.spec,'pumpprobe',false) 
title('Absorptive spectrum')

figure,clf,my2dPlot(s.w,s.w,s.diagrams(1).R,'pumpprobe',false) 
title('Exact R_r')

figure,clf,my2dPlot(s.w,s.w,s.diagrams(2).R,'pumpprobe',false) 
title('Exact R_{nr}')

peak_pos = [2350 2350]; %w1 w3 position
[v,V] = inhomogeneityIndex(s.w,s.w,s.diagrams(1).R,s.diagrams(2).R,peak_pos);


figure,clf,my2dPlot(s.w,s.w,sin(pi/2*V),'pumpprobe',false) 
title('Xtracted c(t)/c(0)')

%% Look at R_r and R_nr in time domain
ss = s;
n_t = s.n_t;
tt = 0:(n_t-1);
ss = ss.calcDiagramsTime;
figure(4),clf,my2dPlot(tt,tt,ss.diagrams(1).R,'pumpprobe',false) 
title('Exact R_r(t1,t2,t3)')

figure(5),clf,my2dPlot(tt,tt,ss.diagrams(2).R,'pumpprobe',false) 
title('Exact R_{nr}(t1,t2,t3)')


%% try to isolate similar parts in ifft of abs spectrum
spec = s.spec;
spec = ifftshift(spec);

n_t = aRFoptions.n_zp/2;
t = 0:(2*n_t-1);
w = s.w;

S1 = ifft2(spec);
figure(101),clf,my2dPlot(t,t,real(S1),'pumpprobe',false) 
%%

S1 = fliplr(circshift(S1,[0 -1]));
S2 = ifft2(spec);

S1(n_t+1:end,:) = 0;
S1(:,n_t+1:end) = 0;

S2(:,n_t+1:end) = 0;
S2(n_t+1:end,:) = 0;

figure(101),clf,my2dPlot(t,t,real(S1),'pumpprobe',false) 
figure(102),clf,my2dPlot(t,t,real(S2),'pumpprobe',false) 
%%
S1 = sgrsfft2(S1);
S1 = fliplr(circshift(S1,[0 -1]));
S1 = fftshift(S1);

S2 = sgrsfft2(S2);
S2 = fftshift(S2);

figure(103),clf,my2dPlot(w,w,real(S1),'pumpprobe',false) 
title('Extracted R_r')

figure(104),clf,my2dPlot(w,w,real(S2),'pumpprobe',false) 
title('Extracted R_{nr}')


%% calculate R and NR from absorptive spectrum
% defaults
options.n_w = 128; %number of freq points
options.phase = 0;
options.flag_plot = true;
%opt.range1 = [w3(1) w3(end)];
%opt.range3 = [w3(1) w3(end)];
options.range1 = [2250, 2450];
options.range3 = [2250, 2450];

result = fromAbsorptiveToRandNR(s.w,s.w,s.spec,options);

%% phasing
options.peak_pos =[2350, 2350];
options.flag_plot = true;
options.range1 = [2250, 2450];
options.range3 = [2250, 2450];

phase = intrinsicPhasing(result.w1,result.w3,result.S,options)

%% inhomogeneity index
peak_pos = [2350, 2350]; %w1 w3 position
[v,V] = inhomogeneityIndex(result.w1,result.w3,result.R,result.NR,peak_pos);
figure,clf,my2dPlot(result.w1,result.w3,sin(pi/2*V),'pumpprobe',false) 
title('Extracted c(t)/c(0)')

%% weakly anharmonic oscillator
clear aRFoptions s lsf lsfParams dyn dynOther
lsfParams = fitParamBnd('Delta1_cm',20,10,30,'');
lsfParams(2) = fitParamBnd('tau1',2,0.5,5,'');
lsf = lsf1expBnd(lsfParams,'free');

on = @(T1,t2,T3,p) 1;
off = @(T1,t2,T3,p) 0;

label = 'ones and zeros'         
fun_array = {on, off};
ind_array = {[1:6], [7:8]};
dyn = additionalDynamics(label,fun_array,ind_array);
aRFoptions.dynOther = dyn; 

aRFoptions.dt = 0.01;
aRFoptions.n_t = 256;
aRFoptions.n_zp = 2*aRFoptions.n_t;
aRFoptions.flag_rotating_frame = true;
aRFoptions.mu01sq = fitParamBnd('mu01sq',1,.1,2,'fixed');
aRFoptions.w_01_cm = fitParamBnd('w_01_cm',2350,2340,2360,'fixed');
aRFoptions.anh_cm = fitParamBnd('anh_cm',50,49,51,'fixed');
aRFoptions.mu12sqRatio = fitParamBnd('mu12sqRatio',2,1.9,2.1,'fixed');
aRFoptions.phase_deg = fitParamBnd('phase_deg',0,-10,20,'fixed');
aRFoptions.damping = lsf;
aRFoptions.t2_array = 0;
aRFoptions.w1_in = 2200:10:2400;
aRFoptions.w3_in = 2200:10:2400;
        


s = aRFWAOBnd(aRFoptions);
s = s.calcSpectrum(s.p0);

figure,clf,my2dPlot(s.w,s.w,s.spec,'pumpprobe',false) 
title('Absorptive spectrum')

%% calculate R and NR from absorptive spectrum
% defaults
options.n_w = 128; %number of freq points
options.phase = 0;
options.flag_plot = true;
%opt.range1 = [w3(1) w3(end)];
%opt.range3 = [w3(1) w3(end)];
options.range1 = [2250, 2450];
options.range3 = [2250, 2450];

result = fromAbsorptiveToRandNR(s.w,s.w,s.spec,options);

%% inhomogeneity index
peak_pos = [2350, 2350]; %w1 w3 position
[v,V] = inhomogeneityIndex(result.w1,result.w3,result.R,result.NR,peak_pos);
figure,clf,my2dPlot(result.w1,result.w3,sin(pi/2*V),'pumpprobe',false) 
title('Extracted c(t)/c(0)')

%% phasing
options.peak_pos =[2350, 2350];
options.flag_plot = true;
options.range1 = [2250, 2450];
options.range3 = [2250, 2450];

phase = intrinsicPhasing(result.w1,result.w3,result.S,options)

%% TRY WITH EXPERIMENTAL DATA

%cwd = pwd;
%datadir = '~/OneDrive - University of Pittsburgh/data/2dir_data/2017/2017-10-18';
%cd(datadir)

m = load('003.mat'); %2017-10-18-003
d = m.data(1);
%cd(cwd)

range1 = [2320, 2370];
range3 = [2250, 2450];
ind1 = find(d.w1>range1(1) & d.w1<range1(2));
ind3 = find(d.w3>range3(1) & d.w3<range3(2));

figure(1),clf
my2dPlot(d.w1(ind1),d.w3(ind3),d.R(ind3,ind1),'pumpprobe',false)
title(sprintf('t_2 = %8.1f ps',d.t2/1000))


% defaults
options.n_w = 256; %number of freq points
options.phase = 0;
options.flag_plot = true;
%opt.range1 = [w3(1) w3(end)];
%opt.range3 = [w3(1) w3(end)];
options.range1 = range1;
options.range3 = range3;

result = fromAbsorptiveToRandNR(d.w1,d.w3,d.R,options);

%% try phasing

options.peak_pos =[2340.4, 2335.1; 2340.4, 2308.8];
options.flag_plot = true;
options.flag_plot = true;
options.range1 = [2320, 2370];
options.range3 = [2250, 2400];

phase = intrinsicPhasing(result.w1,result.w3,result.S,options)

%% inhomogeneity index
% try with function
peak_pos = [2.340830964899873, 2.335883365565500]*1e3; %w1 w3 position
[v,V] = inhomogeneityIndex(result.w1,result.w3,result.R,result.NR,peak_pos);
figure,clf,my2dPlot(result.w1,result.w3,sin(pi/2*V),'pumpprobe',false) 
title('Extracted c(t)/c(0)')
