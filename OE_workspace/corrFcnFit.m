function [output] = corrFcnFit(t2_array,corrStruct,CFoptions)
% I'm assuming the inputs need to be a list of t2 times and a structure
% that contains the correlation values we've previously been extracting.
% Maybe also some sort of functional form, and bounds on the fit.

ub = CFoptions.ub;
lb = CFoptions.lb;
startpoint = CFoptions.startpoint;
corrFcn = CFoptions.fitfcn;
flag_plot = CFoptions.flag_plot;

% extract correlation values and confidence intervals
c2 = zeros(size(t2_array));
c2_std = c2;
for ii = 1:length(t2_array)
    c2(ii) = corrStruct(ii).fitresult.m;
    dummy = confint(corrStruct(ii).fitresult);
    c2_std(ii) = (dummy(2,1) - dummy(1,1))/2;
end

wt = 1./(c2_std).^2;
wt_mean = mean(wt);
wt = wt./wt_mean;

% choose a functional form

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
