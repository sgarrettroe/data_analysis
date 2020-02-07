function [output] = corrFcnFit(t2_array,c2,c2_std,CFoptions)
% CFFitStructure = corrFcnFit(t2_array,c2,c2_std,CFoptions);
%
% This function will fit a set of population times (t2 array) to 
% correlations, weighting the individual points by the inverse of their
% variance (1/std^2). It is currently set up to do a robust (Bisquare) fit.
% Your fitting function can be of any form, though in most cases it should
% be an/(a set of) exponential decay(s).
%
% CFoptions should be a structure with the following fields:
%       # A fitting function ('fitfcn')
%       # A starting point for the fit ('startpoint')
%       # A lower bound for the fit ('lb')
%       # An upper bound for the fit ('ub')
%       # An optional 'flag_plot' field, which, if set to 1, will plot the
%         resulting fit against the data points (as errorbars).
%
% TODO: Allow a variable-length input argument list to define whether you
% want a robust fit, and if so, what type.


ub = CFoptions.ub;
lb = CFoptions.lb;
startpoint = CFoptions.startpoint;
corrFcn = CFoptions.fitfcn;

if isfield(CFoptions,'flag_plot')
    flag_plot = CFoptions.flag_plot;
else
    flag_plot = 0;
end

wt = 1./(c2_std).^2;
wt_mean = mean(wt);
wt = wt./wt_mean;

% do the fitting
[fitresult,gof,fitinfo] = fit(t2_array(:), c2(:), corrFcn, 'StartPoint', startpoint,...
    'lower', lb,'upper',ub,'Weight',wt,'Robust','Bisquare');

output.fitresult = fitresult;
output.gof = gof;
output.fitinfo = fitinfo;
output.t2 = t2_array;
output.c2 = c2;
output.c2_std = c2_std;

if flag_plot
    figure,clf
    hold on
    errorbar(t2_array,c2,c2_std,'rx')
    plot(t2_array,fitresult(t2_array),'k')
    xlim([-5 max(t2_array)+5])
    hold off
    xlabel('t_2 (ps)')
    ylabel('c_2')
end