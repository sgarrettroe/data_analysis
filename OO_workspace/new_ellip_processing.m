function [ellipticityStruct,CFFit] = new_ellip_processing(data,options,CFoptions)

t2_array = [data.t2]./1000;
range1 = options.range1;
range3 = options.range3;

if isfield(options,'flag_plot')
    flag_plot = options.flag_plot;
else
    flag_plot = 0;
end
CFoptions.flag_plot = flag_plot;
% flag_plot will later determine if we're going to plot our CLS fits
% against the data that generated them.

data = cropData(data,range1,range3);
ellipticityStruct = fitEllip(data,options);
[ellip,std_e] = c2FromEllip(t2_array,ellipticityStruct);
CFFit = corrFcnFit(t2_array,ellip,std_e,CFoptions);