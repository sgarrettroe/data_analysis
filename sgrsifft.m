function x2 = sgrsifft(x,varargin)
%sgrsifft.m Calculate the fft but average the first and last points first 
%to remove the baseline artifacts from having a discontinuous 
%function.
%usage:
% x2 = sgrsifft(x)
%
% does a 1D fft of data with the first point replaced by the average of the
% first and last points
%
% x2 = sgrsifft(x,n_zp)
%
% additionally zero pads both dimensions to n_zp in length
%
% x2 = sgrsifft(x,[],dim) or sgrsfft(x,n_zp,dim) applies the fft operation
% across the dimension dim.

%the length of the input
n_zp = length(x);
flag_dim = false;
dim = 0;

%look for zero-padding, if any
if nargin>=2
  if ~isempty(varargin{1})
    n_zp =varargin{1};
  end
end
%look for a dimension specification
if nargin>=3
  flag_dim = true;
  dim = varargin{2};
end

%copy the input (maybe not necessary but this keeps the steps similar to
%the 2d algorithm, where it is necessary)
x2 = x;

%calculate the average and replace the first point
x2(1) = (x(1)+x(end))/2;

%to the fft with any zero-padding
if flag_dim
  x2 = ifft(x2,n_zp,dim);
else
  x2 = ifft(x2,n_zp);
end
