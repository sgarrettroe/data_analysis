function my3dSlices2(w1,w3,w5,R,ind1,varargin)
%my3dSlices2.m 
% call as my3dSlices2(w1,w3,w5,R,ind1)
%
%see also my3dSlices for displaying the 3d object with cuts indicated
n_contours = 12;
flag_fwhm  = false;
while length(varargin)>=2
  arg = varargin{1};
  val = varargin{2};
  switch lower(arg)
    case 'n_contours'
      n_contours = val;
    case 'fwhm'
      if strcmpi(val,'on')
        flag_fwhm = true;
      else
        flag_fwhm = false;
      end
    otherwise
      error(['my3dSlices2: unknown option ',arg])
  end
  varargin = varargin(3:end);
end

s=squeeze(R(:,ind1,:))';
MAX = max(max(abs(s(:))));
level_list = linspace(-MAX,MAX,n_contours+2);
%myCaxis(s);
dl = level_list(2)-level_list(1);
cmin =level_list(1)-dl; %-eps('single')
cmax =level_list(end);%-dl/2+10*eps('single')

contourf(w3,w5,s,level_list);
caxis([cmin cmax]);
if flag_fwhm
  hold on
  %keyboard
  contour(w3,w5,s,[cmax/2 cmax/2],'LineWidth',2,'color',[0.01 0 0]);
  hold off
end

l=line([w3(1) w3(end)],[w3(1),w3(end)]);
set(l,'Color',[0 0 0]);

map = myMapRGB2(n_contours+1);
%map = myMapRGB2(65);
colormap(map);
%myCaxis(s);
axis equal tight
xlabel('\omega_3 / 2\pic')
ylabel('\omega_5 / 2\pic')

%place text of the value of w1
dx = w3(end)-w3(1);
dy = w5(end)-w5(1);
x0 = w3(1)+0.1*dx;
y0 = w5(end)-0.12*dy;
text(x0,y0,['\omega_1 = ' sprintf('%4.0f',w1(ind1)), ' cm^{-1}']);
set(gca,'XLim',[w3(1) w3(end)],'YLim',[w5(1) w5(end)])
set(gca,'TickDir','out')
