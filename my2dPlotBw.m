function varargout=my2dPlot(x,y,z,n_contours)
%my2dPlot should plot contour plots and projections of 2d data
%returning handles to the individual axes that make the plot

clf
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


contour(x,y,z,n_contours,'k');

%axis square%this makes things not line up well...
a(1)=gca;
set(a(1),'Position',[cont_left cont_bot cont_w cont_h]);
line([x(1) x(end)],[x(1) x(end)],'Color',[0 0 0]);
xlim = get(a(1),'XLim');
ylim = get(a(1),'YLim');

a(2)=axes('Position',[cont_left, cont_bot+cont_h, cont_w, proj_h]);
line(x,sum(z,1),'Color',[0 0 0]);
set(a(2),'XTickLabel',[],'YTickLabel',[],'Xlim',xlim);

a(3)=axes('Position',[cont_left+cont_w, cont_bot, proj_h, cont_h]);
line(sum(z,2),y,'Color',[0 0 0]);
set(a(3),'XTickLabel',[],'YTickLabel',[],'YLim',ylim);

for i = 1:nargout
  varargout{i} = a(i);
end
