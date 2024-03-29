function x2 = sgrsfft2(x,varargin)
%sgrsfft2.m Calculate the two-dimensional fft and account for the heaviside
%step function in a single sided fourier transform by averaging the first
%and last points (row and column) of a matrix.
%
%usage:
% x2 = sgrsfft2(x)
%
% does a 2D fft of data with the first and last columns averaged
%
% x2 = sgrsfft2(x,n_zp)
%
% additionally zero pads both dimensions to n_zp in length

%check the size of the input 
n_zp1 = size(x,2);
n_zp3 = size(x,1);
% look for any zero-padding
if nargin>=2
  n_zp1 = varargin{1};
  n_zp3 = n_zp1; 
end

%initialize the intermediate
x2 = x;
%x2(:,1) = (x(:,1)+x(:,end))./2;
%x2(1,2:end) = (x(1,2:end)+x(end,2:end))./2;
%x2 = fft2(x2,n_zp,n_zp);

%do the averaging of the first dimension
x2(1,:) = (x(1,:)+x(end,:))./2;

%now fft that dimension with zero padding if any
x2 = fft(x2,n_zp3);

%now average the second dimension (in frequency space)
x2(:,1) = (x2(:,1)+x2(:,end))./2;

%the second dimension fft
x2 = fft(x2.',n_zp1).';
