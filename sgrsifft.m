function x2 = sgrsifft(x,varargin)
%sgrsifft.m Calculate the fft but average the first and last points first 
%to remove the baseline artifacts from having a discontinuous 
%function.
%usage:
% x2 = sgrsfft(x)
%
% does a 1D fft of data with the first point replaced by the average of the
% first and last points
%
% x2 = sgrsfft(x,n_zp)
%
% additionally zero pads both dimensions to n_zp in length

%the length of the input
n_zp = length(x);

%look for zero-padding, if any
if nargin==2
  n_zp =varargin{1};
end

%copy the input (maybe not necessary but this keeps the steps similar to
%the 2d algorithm, where it is necessary)
x2 = x;

%calculate the average and replace the first point
x2(1) = (x(1)+x(end))/2;

%to the fft with any zero-padding
x2 = ifft(x2,n_zp);
