function varargout = rb2dPlot(varargin)
% rb2dPlot
%
% 20110528, RB: started the function
%
% INPUT:
% either: a struct, it will plot the 2D spectrum
% or: x_axis, y_axis, data, it will plot this
% the x_axis and y_axis should correspond to the size of the matrix
% CONTRARY TO my2dPlot THESE ARE THE WHOLE AXIS, NOT THE PART YOU WANT TO
% PLOT!!! 
%
% other INPUT: 
% - xlim, ylim (opt): the limits for the x and y axes. [0 0] means
% the whole range will be plotted, [0 -1] means that the other axis will be
% used for the range and [n m] means that range will be plotted. In case of
% errors (out of range etc) it will fall back to plot everything. Default
% is xlim = [0 0] and ylim = [0 -1]. 
% - zlimit (opt), number: will change the intensity range. 0 means
% everything is plotted (default). 0 < zlimit <= 1 will plot a range.
% zlimit > 1 will plot an absolute value. The zlimit is always symmetric
% around 0. 
% - n_contours (opt), number: the amount of contours. Using an odd number will give a
% warning. Default is 12.
% - pumprobe (opt), BOOL: changes the axes. Default is FALSE.
% - title (opt), string: a title will be given to the plot.
%



% set variables
n_contours = 12;
zlimit = 0;
flag_pumpprobe = false;
xlim = [0 0];
ylim = [0 -1];
title_string = '';
flag_no_units = false;
x_label = '\omega_3 / 2\pic';
y_label = '\omega_1 / 2\pic';

disp(varargin)

% read varargin (part I)
if isa(varargin{1}, 'struct')
  x_axis = varargin{1}.w3;
  y_axis = varargin{1}.w1;
  data = varargin{1}.R;
  varargin = varargin(2:end);
else
  x_axis = varargin{1};
  y_axis = varargin{2};
  data = varargin{3};
  varargin = varargin(4:end);
end

% read varargin
while length(varargin) >= 2
  arg = varargin{1};
  val = varargin{2};
  
  switch lower(arg)
    case 'xlim'
      xlim = val;
    case 'ylim'
      ylim = val;
    case 'zlimit'
      zlimit = val;
    case 'n_contours' 
      n_contours = val;
    case 'pumpprobe'
      flag_pumpprobe = val;
    case 'title'
      title_string = val;
    case 'no_units'
      flag_no_units = val;
    otherwise
      error(['rb2dPlot: unknown option ', arg]);
  end 
  varargin = varargin(3:end);
end

% error checking
if mod(n_contours,2)
  warning('my2dPlot4: Odd number of contour lines may produce unexpected results!')
end

% determine the x and y axes
if xlim == [0 0]
  xrange = find(x_axis);
  if ylim == [0 0]
    yrange = find(y_axis);
  elseif ylim == [0 -1]
    yrange = find(y_axis >= x_axis(1) & y_axis <= x_axis(end));
  else 
    yrange = find(y_axis >= ylim(1) & y_axis <= ylim(end));
  end
elseif xlim == [0 -1]
  if ylim == [0 0]
    yrange = find(y_axis);
    xrange = find(x_axis >= y_axis(1) & x_axis <= y_axis(end));
  elseif ylim == [0 -1]
    yrange = find(y_axis);
    xrange = find(x_axis);
  else 
    yrange = find(y_axis >= ylim(1) & y_axis <= ylim(end));
    xrange = find(x_axis >= y_axis(1) & x_axis <= y_axis(end));
  end 
else
  xrange = find(x_axis >= xlim(1) & x_axis <= xlim(end));
  if ylim == [0 0]
    yrange = find(y_axis);
  elseif ylim == [0 -1]
    yrange = find(y_axis >= xlim(1) & y_axis <= xlim(end));
  else 
    yrange = find(y_axis >= ylim(1) & y_axis <= ylim(end));
  end 
end
  
% set the x, y, and z
if flag_no_units
  x = xrange;
  y = yrange;
else
  x = x_axis(xrange);
  y = y_axis(yrange);
end
z = data(yrange,xrange);

% load the color scheme
map = myMapRGB2(n_contours);

% determine the range to be plotted
if zlimit <= 0 
  % 0: use the whole range
  [ca, level_list]= myCaxis2(z, n_contours);
elseif zlimit > 0 && zlimit <= 1
  % 0 < zlimit <= 1: plot a ratio
  [ca, level_list] = myCaxis2(z, n_contours);
  ca = ca * zlimit;
  level_list = level_list * zlimit;
else
  % otherwise, plot between abs. values
  ca = [-zlimit zlimit];
  level_list = linspace(-zlimit, zlimit, n_contours+2);
end

disp(ca);

% plot, use colormap and set axes
contourf(x, y, z, level_list);
colormap(map);
caxis(ca);

% labels for the plot
if flag_pumpprobe
  x_label = '\omega_{probe} / 2\pic';
  y_label = '\omega_{pump} / 2\pic';
elseif flag_no_units
  x_label = 'pixels';
  y_label = 'FT(time)';
%else
%  x_label = '\omega_3 / 2\pic';
%  y_label = '\omega_1 / 2\pic';
end
xlabel(x_label);
ylabel(y_label);
title(title_string);

%a(1) = gca;
%set(a(1), 'Position', 



%plot([0 0], [3000 3000])
%line([0 0], [3000 3000])


