function varargout = rb2dPlot(s, varargin)

% set variables
n_contours = 12;
zlimit = 0;
flag_pumpprobe = false;
xlim = [0 0];
ylim = [0 -1];
title_string = '';

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
    otherwise
      error(['rb2dPlot: unknown option ', arg]);
  end 
  varargin = varargin(3:end);
end

% determine the x and y axes
if xlim == [0 0]
  xrange = find(s.w3);
  if ylim == [0 0]
    yrange = find(s.w1);
  elseif ylim == [0 -1]
    yrange = find(s.w1 >= s.w3(1) & s.w1 <= s.w3(end));
  else 
    yrange = find(s.w1 >= ylim(1) & s.w1 <= ylim(end));
  end
elseif xlim == [0 -1]
  if ylim == [0 0]
    yrange = find(s.w1);
    xrange = find(s.w3 >= s.w1(1) & s.w3 <= s.w1(end));
  elseif ylim == [0 -1]
    yrange = find(s.w1);
    xrange = find(s.w3);
  else 
    yrange = find(s.w1 >= ylim(1) & s.w1 <= ylim(end));
    xrange = find(s.w3 >= s.w1(1) & s.w3 <= s.w1(end));
  end 
else
  xrange = find(s.w3 >= xlim(1) & s.w3 <= xlim(end));
  if ylim == [0 0]
    yrange = find(s.w1);
  elseif ylim == [0 -1]
    yrange = find(s.w1 >= xlim(1) & s.w1 <= xlim(end));
  else 
    yrange = find(s.w1 >= ylim(1) & s.w1 <= ylim(end));
  end 
end
  
% set the x, y, and z
x = s.w3(xrange);
y = s.w1(yrange);
z = s.R(yrange,xrange);

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

if flag_pumpprobe
  x_label = '\omega_{probe} / 2\pic';
  y_label = '\omega_{pump} / 2\pic';
else
  x_label = '\omega_1 / 2\pic';
  y_label = '\omega_3 / 2\pic';
end
xlabel(x_label);
ylabel(y_label);
title(title_string);



