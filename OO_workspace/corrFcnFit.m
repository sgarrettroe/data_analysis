function [output] = corrFcnFit(t2_array,c2,c2_std,CFoptions)
% I'm assuming the inputs need to be a list of t2 times and a structure
% that contains the correlation values we've previously been extracting.
% Maybe also some sort of functional form, and bounds on the fit.

ub = CFoptions.ub;
lb = CFoptions.lb;
startpoint = CFoptions.startpoint;
corrFcn = CFoptions.fitfcn;
flag_plot = CFoptions.flag_plot;

wt = 1./(c2_std).^2;
wt_mean = mean(wt);
wt = wt./wt_mean;

% do the fitting
[fitresult,gof,fitinfo] = fit(t2_array(:), c2(:), corrFcn, 'StartPoint', startpoint,...
    'lower', lb,'upper',ub,'Weight',wt,'Robust','Bisquare');

output.fitresult = fitresult;
output.gof = gof;
output.fitinfo = fitinfo;

if flag_plot
    figure(1),clf
    hold on
    errorbar(t2_array,c2,c2_std,'rx')
    plot(t2_array,fitresult(t2_array))
    hold off
    xlabel('t_2 (ps)')
    ylabel('CLS (c_2)')
end