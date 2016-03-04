function [ellip,std_e] = c2FromEllip(t2_array,ellipticityStruct)
% this will fall apart if we used a different fitting functional form

ellip = zeros(size(t2_array));
std_e = zeros(size(t2_array));
for ii = 1:length(ellipticityStruct)
  dummy = confint(ellipticityStruct(ii).fitresult);
  err_sD = dummy(:,4);
  err_sA = dummy(:,5);
  std_sD = (err_sD(2) - err_sD(1))/2; %std of sd
  std_sA = (err_sA(2) - err_sA(1))/2; %std of sa
  sD = ellipticityStruct(ii).fitresult.sD;
  sA = ellipticityStruct(ii).fitresult.sA;
  de_sD = (4*sD*sA^2)/(sD^2+sA^2)^2; %derivatives for sD
  de_sA = -(4*sD^2*sA)/(sD^2+sA^2)^2; %derivatives for sA
  ellip(ii) = (sD.^2-sA.^2)./(sD.^2+sA.^2);
  std_e(ii) = (de_sD.^2*std_sD.^2 + de_sA.^2*std_sA.^2)^0.5;
end
