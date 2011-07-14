function varargout = myCaxis2(Z, n_contours)
% myCaxis will determine the levels for the 2d plot. 
%
% INPUT
% - Z (2d array): the data
% - n_contours (number): the amount of contours needed. When using an even
%   number, the data will be centered around 0. With an odd number, you
%   might run into problems. This will give an error.
%
% OUTPUT
% - [cmin cmax]: the upper and lower end of the level list
% - level_list (array): an array with n_contours+1 numbers. 
%
% NOTE
% I'm not sure why myCaxis2 was so complicated.

if mod(n_contours,2) == 0
  warning('myCaxis2: Odd number of contour lines may produce unexpected results!')
end

MAX = max(abs(Z(:)));

level_list = linspace(-MAX, MAX, n_contours + 1);

cmin = level_list(1);
cmax = level_list(end);

if nargout>=1
  varargout{1} = [cmin cmax];
end
if nargout>=2
  varargout{2} = level_list;
end
