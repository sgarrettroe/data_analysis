function x3 = sgrsifft3(x,varargin)
%this does the 3d inverse FFT on a tensor and reduces the artifact which comes from
%the a single sided nature of the data (i.e. it isn't periodic). So the
%function first averages the t=0 planes of the 3D data first _then_ does
%the FFT.
%sided 
n_zp = size(x,1);
if nargin>=2
  n_zp = varargin{1};
end
%x(:,:,1) = (x(:,:,1)-x(:,:,end))./2;
%x(:,1,:) = (x(:,1,:)-x(:,end,:))./2;
%x(1,:,:) = (x(1,:,:)-x(end,:,:))./2;
x(:,:,1) = (x(:,:,1))./2;
x(:,1,:) = (x(:,1,:))./2;
x(1,:,:) = (x(1,:,:))./2;
x3 = ifftn(x,[n_zp,n_zp,n_zp]);
