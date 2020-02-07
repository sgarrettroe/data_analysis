function varargout=my2dPlot(x,y,z,varargin)

%my2dPlot should plot contour plots and projections of 2d data
%returning handles to the individual axes that make the plot
%function varargout=my2dPlot(x,y,z,['opt',val])
%
%RB, 03052011: introduced variable zlimit, which changes the range over 
%which the contours are drawn. zlimit <= 0 is default and uses the data, 
%zlimit > 0 and <= 1 will plot a ratio, zlimit > 0 limits the contours to 
%the range specified.
%
%fix contour levels and colors

n_contours = 12;
flag_pumpprobe = true;
zlimit = 0;
while length(varargin)>=2
  arg = varargin{1};
  val = varargin{2};
  switch lower(arg)
    case 'zlimit'
      zlimit = val;  
    case 'n_contours'
      n_contours = val;
      if mod(n_contours,2)
        warning('my2dPlot4: Odd number of contour lines may produce unexpected results!')
      end
    case {'pumpprobe_style','pumpprobe'}
      flag_pumpprobe = val;
    otherwise
      error(['my2dPlot: unknown option ',arg])
  end
  varargin = varargin(3:end);
end

pos = get(gcf,'Position');
if get(gcf,'WindowStyle')~='docked',
  set(gcf,'Position',[pos(1) pos(2) 415 416],'Color',[0.9 0.9 0.9]);
end
%contour properties
cont_left = 0.18;
cont_bot = 0.15;
cont_w = 0.6;
cont_h = 0.6;

%projecttion properties
proj_h= 0.2;

map=myMapRGB2(n_contours);
%MAX = max([max(z(:)) abs(min(z(:)))]);
%MAX = max(abs(z(:)));
%level_list = linspace(-MAX,MAX,n_contours+2);

% zlimit determines what is what the range is that is covered by the
% contours. If limit <= 0, the data is used to determine it (old behaviour
% and default). If zlimit is between 0 and 1, it will use the data as well,
% but only plots a fraction of it (the number between 0 and 1. If zlimit is
% large than 1, it will use that value.
if zlimit <= 0 
  [ca, level_list]= myCaxis2(z, n_contours);
  %ca 
elseif zlimit > 0 && zlimit <= 1
  [ca, level_list] = myCaxis2(z, n_contours);
  ca = ca * zlimit;
  level_list = level_list * zlimit;
else
  ca = [-zlimit zlimit];
  level_list = linspace(-zlimit, zlimit, n_contours+2);
end
  
  
contourf(x,y,z,level_list);
colormap(map)
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

%axis square%this makes things not line up well...
a(1)=gca;
set(a(1),'Position',[cont_left cont_bot cont_w cont_h]);
line([x(1) x(end)],[x(1) x(end)],'Color',[0 0 0]);
xlim = get(a(1),'XLim');
ylim = get(a(1),'YLim');

a(2)=axes('Position',[cont_left, cont_bot+cont_h, cont_w, proj_h]);
%line(x,sum(z,2)','Color',[0 0 0]);
line(x,sum(z,1),'Color',[0 0 0]);
set(a(2),'XTickLabel',[],'YTickLabel',[],'Xlim',xlim);

a(3)=axes('Position',[cont_left+cont_w, cont_bot, proj_h, cont_h]);
line(sum(z,2)',y,'Color',[0 0 0]);
set(a(3),'XTickLabel',[],'YTickLabel',[],'YLim',ylim);

%linkaxes([a(2) a(1)],'x')
%linkaxes([a(3) a(1)],'y')

hlink_x = linkprop([a(1),a(2)],'XLim');
hlink_y = linkprop([a(1),a(3)],'YLim');
key_x = 'graphics_linkprop';
key_y = 'graphics_linkprop';

setappdata(a(2),key_x,hlink_x);
setappdata(a(3),key_y,hlink_y);
for i = 1:nargout
  varargout{i} = a(i);
end
