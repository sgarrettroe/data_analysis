data = load2DIRdata('c:\Users\Thomas Brinzer\Documents\LocalData\2015-05-19',[9:20 23:29]);
data = sort2DIRdata(data);
%%
view2DIRdata(data,[1950 2120],[1950 2120])
%%
PFOptions.range1 = [2043 2067];
PFOptions.range3 = [1990 2100];
PFOptions.fitfcn = fittype('- a1.*voigt(w,center,wg,wl) + a2.*voigt(w,center - anh,wg,wl)', ...
    'coeff',{'center','a1','a2','wg','wl','anh'},'indep',{'w'});
PFOptions.startpoint = [2055 1000 1000 10 5 24];
PFOptions.lb = [2010 4000 4000 0 0 20];
PFOptions.ub = [2085 Inf Inf 20 5 30];
% PFOptions.zp_factor = 2;
PFOptions.flag_plot = 1;
PFOptions.estimateArea = 1;

CFOptions.startpoint = [0.3 0.3 3 100];
CFOptions.lb = [0 0 0 10];
CFOptions.ub = [1 1 50 500];
CFOptions.fitfcn = fittype(@(a1,a2,t1,t2,x) ...
    a1.*exp(-x./t1) +  a2*exp(-x./t2), ...
    'coeff', {'a1', 'a2', 't1','t2'}, 'indep', {'x'});

[peakFit,CLS,CFFit1] = new_CLS_processing(data,PFOptions,CFOptions);
%%
ellipOptions.range1 = [1950 2120];
ellipOptions.range3 = [1950 2120];
% fitting model F2DS
ellipOptions.startpoint = [2050 2050 10 10 30 0]; % this array will be concatenated with -max(max(xx)) 
ellipOptions.upper_bound = [1000 2100  2100 20 20 60  1e2];
ellipOptions.lower_bound = [-1e6 1900  1900 1  1 0 -1e2];
ellipOptions.flag_plot = 1;

CFOptions.startpoint = [0.3 0.3 3 100];
CFOptions.lb = [0 0 0 10];
CFOptions.ub = [1 1 50 500];
CFOptions.fitfcn = fittype(@(a1,a2,t1,t2,x) ...
    a1.*exp(-x./t1) +  a2*exp(-x./t2), ...
    'coeff', {'a1', 'a2', 't1','t2'}, 'indep', {'x'});

[ellip,CFFit2] = new_ellip_processing(data,ellipOptions,CFOptions);

%%
% format: YYYY-MM-DD-<experiment>_summary.mat
save('2015-05-19-01_summary.mat')