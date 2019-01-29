function FWHM_cm = f_Lor(T2_ps)
% f_Lor Full width at half maximum of a Lorentzian.
% FWHM = f_Lor(T2) is the full width at half maximum of a Lorentzian, given
% the dephasing time (T2) in ps.
% 
% Lorentzian function:
%    L(x) = (1/pi)* Gamma/(x^2 + Gamma^2)
%    Gamma = 1/T2
%
% In ps^-1, the FWHM is trivial to calculate:
%    FWHM = 2*Gamma = 2/T2_ps
% 
% To convert from rad*ps^-1 (angular frequency) to cm^-1, we divide by
% 2*pi*c (in appropriate units)
%     
%    FWHM_cm = (pi*c*T2_ps)^(-1)

c = 2.9979E10; % cm/s
wavenumbersToInvPs = c*1E-12;
invPsToWavenumbers=1/wavenumbersToInvPs;
FWHM_cm = invPsToWavenumbers*1./(pi*T2_ps);