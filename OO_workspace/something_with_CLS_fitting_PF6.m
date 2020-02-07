%Peak Fitting Options:
PFOptions.range1 = [2341 2346];
PFOptions.range3 = [2312 2355];
PFOptions.fitfcn = fittype('- a1.*voigt(w,center,wg,wl) + a2.*voigt(w,center - anh,wg,wl)', ...
    'coeff',{'center','a1','a2','wg','wl','anh'},'indep',{'w'});
PFOptions.startpoint = [2343 2000 2000 3 3 24];
PFOptions.lb = [2335 2000 2000 1 1 22];
PFOptions.ub = [2355 Inf Inf 8 8 30];
PFOptions.zp_factor = 2;
PFOptions.flag_plot = 1;
% 
% Correlation Function Options
CFOptions.startpoint = [0.3 0.3 3 100];
CFOptions.lb = [0 0 0 30];
CFOptions.ub = [1 1 50 500];
CFOptions.fitfcn = fittype(@(a1,a2,t1,t2,x) ...
    a1.*exp(-x./t1) +  a2*exp(-x./t2), ...
    'coeff', {'a1', 'a2', 't1','t2'}, 'indep', {'x'});

% CFOptions.startpoint = [0.3 0.3 3 100 0.05];
% CFOptions.lb = [0 0 0 30 0];
% CFOptions.ub = [1 1 50 500 0.3];
% CFOptions.fitfcn = fittype(@(a1,a2,t1,t2,c,x) ...
%     a1.*exp(-x./t1) +  a2*exp(-x./t2)+  c, ...
%     'coeff', {'a1', 'a2', 't1','t2', 'c'}, 'indep', {'x'});

% CFOptions.startpoint = [0.3 0.3 0.3 3  100 500 ];
% CFOptions.lb =         [0   0   0   0  30  0];
% CFOptions.ub =         [1.5 1.5 1.5 50 300 5000];
% CFOptions.fitfcn = fittype(@(a1,a2,a3,t1,t2,t3,x) ...
%     a1.*exp(-x./t1) +  a2*exp(-x./t2)+  a3*exp(-x./t3), ...
%     'coeff', {'a1', 'a2', 'a3','t1','t2', 't3'}, 'indep', {'x'});


[peakFit,CLS,CFFit] = new_CLS_processing(data,PFOptions,CFOptions);