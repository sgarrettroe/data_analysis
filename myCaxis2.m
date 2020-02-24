function varargout = myCaxis2(Z,n_contours)
%function c = myCaxis(Z,n_contours) 
% where Z is the data and n_contours is the number of contour levels!
%this version is trying to get the stupid coloring to work out right

  %take the larger of the max an min values of Z
  %MAX = max( [max(max(Z)) abs( min(min(Z)) )]);
  %MAX = max([max(Z(:)) abs(min(Z(:)))]);
  MAX = max(abs(Z(:)));

  level_list = linspace(-MAX,MAX,n_contours+2);
  dl = level_list(2)-level_list(1); % I don't really understand why I have to do this shift
  level_list(1) = level_list(1);%-eps('single');
  level_list(end) = level_list(end);%+eps('single');  
  cmin =level_list(1)-dl/2;
  cmax =level_list(end);%-dl/2+10*eps('single');

  caxis([cmin cmax]);
  
if nargout>=1
  varargout{1} = [cmin cmax];
end
if nargout>=2
  varargout{2} = level_list;
end

% 
% 
% 
% function varargout = myCaxis2(Z, n_contours)
% % myCaxis will determine the levels for the 2d plot. 
% %
% % INPUT
% % - Z (2d array): the data
% % - n_contours (number): the amount of contours needed. When using an even
% %   number, the data will be centered around 0. With an odd number, you
% %   might run into problems. This will give an error.
% %
% % OUTPUT
% % - [cmin cmax]: the upper and lower end of the level list
% % - level_list (array): an array with n_contours+1 numbers. 
% %
% % NOTE
% % I'm not sure why myCaxis2 was so complicated.
% 
% if mod(n_contours,2)
%   warning('myCaxis2: Odd number of contour lines may produce unexpected results!')
% end
% 
% MAX = max(abs(Z(:)));
% 
% % When viewing the colorbar in a plot, the result seems to be incorrect,
% % but it is actually a problem with the colorbar, not this script (see test
% % script below the code).
% level_list = linspace(-MAX, MAX, n_contours + 2);
% 
% cmin = level_list(1);
% cmax = level_list(end);
% 
% if nargout>=1
%   varargout{1} = [cmin cmax];
% end
% if nargout>=2
%   varargout{2} = level_list;
% end
% 
% 
% 
% 
% % 
% % % This script can be used to test if the colorbar is correct
% % % produce some data
% % x = linspace(-1, 1, 100);
% % y = linspace(-1, 1, 100);
% % 
% % x_axis = linspace(-1, 1, 100);
% % y_axis = linspace(-1, 1, 100);
% % 
% % [X, Y] = meshgrid(x, y);
% % 
% % data = X + Y;
% % 
% % % some variables
% % n_contours = 6; % number of contours
% % brightness = 0.2; % brightness
% % 
% % % make colormap
% % map = myMapRGB2(n_contours, brightness)
% % 
% % % make contourlevels 
% % [ca, level_list] = myCaxis2(data, n_contours)
% % 
% % % plot data 
% % figure(4),clf
% % contourf(x_axis, y_axis, data, level_list)
% % colormap(map)
% % colorbar()