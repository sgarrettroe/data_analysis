function [ellipticity,CFFit] = new_ellip_processing(data,options,CFoptions)

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

ellipticity = fitEllip(data,options);
% this will fall apart if we used a different fitting functional form

ellip = zeros(size(t2_array));
std_e = zeros(size(t2_array));
for ii = 1:length(ellipticity)
  dummy = confint(ellipticity(ii).fitresult);
  err_sD = dummy(:,4);
  err_sA = dummy(:,5);
  std_sD = (err_sD(2) - err_sD(1))/2; %std of sd
  std_sA = (err_sA(2) - err_sA(1))/2; %std of sa
  sD = ellipticity(ii).fitresult.sD;
  sA = ellipticity(ii).fitresult.sA;
  de_sD = (4*sD*sA^2)/(sD^2+sA^2)^2; %derivatives for sD
  de_sA = -(4*sD^2*sA)/(sD^2+sA^2)^2; %derivatives for sA
  ellip(ii) = (sD.^2-sA.^2)./(sD.^2+sA.^2);
  std_e(ii) = (de_sD.^2*std_sD.^2 + de_sA.^2*std_sA.^2)^0.5;
end

CFFit = corrFcnFit(t2_array,ellip,std_e,CFoptions);