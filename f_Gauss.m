function FWHM = f_Gauss(sigma)
% f_Gauss Full width at half maximum of a Gaussian.
% FWHM = f_Gauss(sigma) returns the FWHM of a Gaussian distribution with
% a standard deviation of sigma.
%
% Gaussian function
%     G(x) = a*exp(-(x-x0)^2/(2*sigma^2)
%     FWHM = 2*sigma*sqrt(2 ln(2))
%
% The FWHM will be returned in the same units as the standard deviation
% (e.g. cm-1 --> cm-1)

FWHM = 2*sigma*sqrt(2*log(2));