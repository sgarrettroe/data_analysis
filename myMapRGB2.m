function map = myMapRGB2(n_levels, varargin)
% myMapRGB2(Z) makes a red-white-blue colormap for the 2d-plots.
%
% RB, 20110713: now you can set hwo pale you want to have it! Default
% min_val is now 0.2, this gives a bit brighter result. It used to be 0.5.
%
% There is a problem with the colorbar in Matlab. Use the test script below
% to make sure that the colors and levels are correct.

min_val = 0.2;

while length(varargin) >= 2
  arg = varargin{1};
  val = varargin{2};  
  switch lower(arg)
    case 'brightness'
      min_val = val;
    otherwise
      error(['process2d: unknown option ', arg]);
  end 
  varargin = varargin(3:end);
end

n_2 = floor(n_levels/2);

dx = (1 - min_val)/n_2;
r=[min_val:dx:1, ones(1,n_2)];
g=[min_val:dx:1, 1-dx:-dx:min_val];
b=[ones(1,n_2), 1:-dx:min_val];

map = [r',g',b'];
  
% 
% % This script can be used to test if the colorbar is correct
% % produce some data
% x = linspace(-1, 1, 100);
% y = linspace(-1, 1, 100);
% 
% x_axis = linspace(-1, 1, 100);
% y_axis = linspace(-1, 1, 100);
% 
% [X, Y] = meshgrid(x, y);
% 
% data = X + Y;
% 
% % some variables
% n_contours = 6; % number of contours
% brightness = 0.2; % brightness
% 
% % make colormap
% map = myMapRGB2(n_contours, brightness)
% 
% % make contourlevels 
% [ca, level_list] = myCaxis2(data, n_contours)
% 
% % plot data 
% figure(4),clf
% contourf(x_axis, y_axis, data, level_list)
% colormap(map)
% colorbar()