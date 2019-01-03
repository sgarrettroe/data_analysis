function FWHM = f_Voigt(f_Lor,f_Gauss)
% f_Voigt Full width at half maximum of a Voigt profile
% FWHM = f_Voigt(f_Lor,f_Gauss), where f_Lor and f_Gauss are the Lorentzian
% and Gaussian full widths at half maximum respectively.
%
% J. J. Olivero and R. L. Longbothum. EMPIRICAL FITS TO THE VOIGT LINE
% WIDTH: A BRIEF REVIEW. J. Quant. Spectrosc. Radiat. Transfer, Vol. 17,
% pp. 233-236. Pergamon Press 1977.

FWHM = 0.5346*f_Lor + sqrt(0.2166* f_Lor^2 + f_Gauss^2);