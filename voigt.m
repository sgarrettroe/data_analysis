function [y] = voigt( wavenumberArray,centerLine,widthGauss, widthLorentz )
% voigt  Calculation of the VOIGT profile 
%
%   [y] = voigt( wavenumberArray,centerLine,widthGauss, widthLorentz )
%   The function calculates the Voight profile using the algorithm 
%   kindly provided by Dr. F. Schreier in FORTRAN and rewritten to MATLAB
%   by Dr. N. Cherkasov
% 
%   For more details on algorithm see the publication:
%   F. Schreier: Optimized Implementations of Rational Approximations for the Voigt ane Complex Error Function. 
%   J. Quant. Spectrosc. & Radiat. Transfer, 112(6), 1010–1025, doi 10.1016/j.jqsrt.2010.12.010, 2011. 
%
%
%   INPUT ARGUMENTS
%       wavenumberArray - array 1*N of wavenumbers 
%       centerLine - position of the band center
%       widthGauss - parameter of the width of the Gaussian component (Half-width at half maximum)
%       widthLorentz - parameter of the width of the Lorentzian component (Half-width at half maximum)
%
% 	OUTPUT
%       y - array 1*N of intensities
%
% The function was used for the deconvolution of IR spectra
% see the publication
%
% 27-December-2013 N. Cherkasov
% Comments and questions to: n.b.cherkasov@gmail.com


% converting to dimensionless coordinates
x=sqrt(log(2)).*(wavenumberArray-centerLine)./(widthGauss);
y=sqrt(log(2)).*(widthLorentz/widthGauss);

w=complexErrorFunction(x,y);
y=sqrt(log(2)/pi)/widthGauss.*real(w);

end

