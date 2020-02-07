function varargout = my3dPlot(w1,w3,w5,v,varargin)
%function varargout = my3dPlot(w1,w3,w5,v,['opt',val])
%this one does orthographic or perspective and I'm trying to get the colors right on the
%projections
%
%opt can be:
% 'labels' = [val1 val2 ...] specify axis labels, otherwise the default set
%                      [2300 2500 2700] is used
% 'ticks'
% 'clf' = {'on','off'} determines whether or not to clear the current
%                      figure (turn off for subplots)
% 'xlabel' = {'on','off'}
% 'projection' = {'orthographic','projection'}
% 'n_contours' = val the number of plotted contours
% 'n_shells' = {2,3} the number of shells plotted
% 'contour_plots' = {|'on'|,'off'} do or do not plot the xyz projections
% 'max_method' = {|'abs'|,'pos','neg'} how to normalize the spectrum
% 'alpha' = [0.1 0.5 0.75] vector of 3 transparency values
% 'iso_value' = [0.1 0.25 0.5]; %vector of 3 isosurface values

%pos = get(gcf,'Position')'
%set(gcf,'Position',[pos(1) pos(2) 712 553]);
iso_value = [0.1 0.25 0.5];
alpha_value = [0.1 0.5 0.75];
label_mode = 'auto';
tick_mode = 'auto';
labels=[2350 2500 2650];
ticks= [2350 2500 2650];
n_contours = 12;
n_shells = 2;
projection_type = 'orthographic';
flag_clf = true;
flag_xlabel = true;
flag_ylabel = true;
flag_zlabel = true;
flag_proj = true;
max_method = 'abs';
while length(varargin)>=2
  arg = varargin{1};
  val = varargin{2};
  switch lower(arg)
    case {'labels','label'}
      label_mode = 'man';
      labels = val;
    case 'ticks'
      tick_mode = 'man';
      ticks = val;
    case 'clf'
      if strcmpi(val,'on')
        flag_clf = true;
      elseif strcmpi(val,'off')
        flag_clf = false;
      else
        error(['unknown value for option ''clf'' = ',val]);
      end
    case 'xlabel'
      if strcmpi(val,'on')
        flag_xlabel = true;
      elseif strcmpi(val,'off')
        flag_xlabel = false;
      else
        error(['unknown value for option ''xlabel'' = ',val]);
      end
    case 'ylabel'
      if strcmpi(val,'on')
        flag_ylabel = true;
      elseif strcmpi(val,'off')
        flag_ylabel = false;
      else
        error(['unknown value for option ''ylabel'' = ',val]);
      end
    case 'zlabel'
      if strcmpi(val,'on')
        flag_zlabel = true;
      elseif strcmpi(val,'off')
        flag_zlabel = false;
      else
        error(['unknown value for option ''zlabel'' = ',val]);
      end
    case 'n_contours'
      n_contours = val;
      if mod(n_contours,2)
        warning('my3dPlot3: Odd number of contour lines may produce unexpected results!')
      end
    case 'n_shells'
      if val == 2 || val == 3
        n_shells = val;
      else
        warning(['n_shells can only be 2 or 3, not ' num2str(val)]);
      end
    case 'alpha'
      if length(val)<3,
        error('alpha must be a 3 element vector');
      end
      alpha_value = val;
    case {'contour_plots','contour_plots'}
      if strcmpi(val,'on')
        flag_proj = true;
      elseif strcmpi(val,'off')
        flag_proj = false;
      else
        error(['my3dPlot3: ''contour_plot'' must be on or off not ' val]);
      end
    case {'projection','projection_type'}
      projection_type = val;
    case {'max_method','norm_method'}
      max_method = val;
    case 'iso_value'
      iso_value = val;
    otherwise
      error(['my3dPlot3: unknown option ',arg])
  end
  varargin = varargin(3:end);
end

%clear current figure
if flag_clf,clf,end 

%if they are matrices, make them row vectors
if ndims(w1)==3
    w1 = reshape(squeeze(w1(1,:,1)),1,[]);
end
if ndims(w3)==3
    w3 = reshape(squeeze(w3(:,1,1)),1,[]);
end
if ndims(w5)==3
    w5 = reshape(squeeze(w5(1,1,:)),1,[]);
end

switch max_method
    case 'abs'
      v = v./max(abs(v(:)));
    case 'pos'
      v = v./max(v(:));
    case 'neg'
      v = v./max(-v(:));
  end
  
  %v = v./max(max(max(v)));
  %map = myMapRGB2(n_contours+1);
  map = myMapRGB2(n_contours);
  colormap(map);
  
  hold on

  if flag_proj %if using projections
    %look for maximum in projections to draw contours
    dx = w1(2)-w1(1);
    dy = w3(2)-w3(1);
    dz = w5(2)-w5(1);
    s=squeeze(sum(v,2))'*dx; %yz plane
    MAX = max(abs(s(:)));
    s=squeeze(sum(v,1))'*dy; %xz plane
    tmp = max(abs(s(:)));
    if tmp>MAX, MAX = tmp;end
    s=sum(v,3)*dz; %xy plane
    tmp = max(abs(s(:)));
    if tmp>MAX, MAX = tmp;end

    %now compute level list
    %should use myCaxis2 here...
    level_list = linspace(-MAX,MAX,n_contours+2);
    dl = level_list(2)-level_list(1); % I don't really understand why I have to do this shift
    cmin =level_list(1)-dl/2;
    cmax =level_list(end);%-dl/2+10*eps('single');
    

    
    %yz plane
    s=squeeze(sum(v,2))'*dx;
    [c,h]=contourf(w3,w5,s,level_list);
    if cmax-cmin>0,
      caxis([cmin cmax]);
    end
    %myCaxis(s);
    %colormap(map),myCaxis(s);
    h=findobj(h,'type','patch');
    for j=1:length(h),
      xd=w1(end).*ones(size(get(h(j),'YData')));
      zd=get(h(j),'YData');
      yd=get(h(j),'XData');
      set(h(j),'XData',xd,'YData',yd,'ZData',zd);
    end
    set(h,'FaceLighting','none');

    %xz plane
    s=squeeze(sum(v,1))'*dy;
    [c,h]=contourf(w1,w5,s,level_list);
%    caxis([cmin cmax]);
    %colormap(map),myCaxis(s);
    h=findobj(h,'type','patch');
    for j=1:length(h),
      xd=get(h(j),'XData');
      zd=get(h(j),'YData');
      yd=w3(end).*ones(size(get(h(j),'YData')));
      set(h(j),'XData',xd,'YData',yd,'ZData',zd);
    end
    set(h,'FaceLighting','none');
    
    %xy plane
    s=sum(v,3)*dz;
    [c,h]=contourf(w1,w3,s,level_list);
%    caxis([cmin cmax]);
    h=findobj(h,'type','patch');
    %colormap(map),myCaxis(s);
    for j=1:length(h),zd=w5(1).*ones(size(get(h(j),'XData')));set(h(j),'ZData',zd);end
    set(h,'FaceLighting','none');
    
  end %end of projections
  
  if n_shells>2
    p10 = patch(isosurface(w1,w3,w5,v,iso_value(1)));
  end
  p25 = patch(isosurface(w1,w3,w5,v,iso_value(2)));
  p50 = patch(isosurface(w1,w3,w5,v,iso_value(3)));
  %isonormals(w1,w3,w5,v,p10);
  if n_shells>2
    p10m = patch(isosurface(w1,w3,w5,v,-iso_value(1)));
  end
  p25m = patch(isosurface(w1,w3,w5,v,-iso_value(2)));
  p50m = patch(isosurface(w1,w3,w5,v,-iso_value(3)));
  %isonormals(w1,w3,w5,v,p25);
  isonormals(w1,w3,w5,v,p50);
  if n_shells>2
    set(p10,'FaceColor','red','EdgeColor','none','FaceAlpha',alpha_value(1));
  end
  set(p25,'FaceColor','red','EdgeColor','none','FaceAlpha',alpha_value(2));
  set(p50,'FaceColor','red','EdgeColor','none','FaceAlpha',alpha_value(3));
  if n_shells>2
    set(p10m,'FaceColor','blue','EdgeColor','none','FaceAlpha',alpha_value(1));
  end
  set(p25m,'FaceColor','blue','EdgeColor','none','FaceAlpha',alpha_value(2));
  set(p50m,'FaceColor','blue','EdgeColor','none','FaceAlpha',alpha_value(3));
  ha=gca; 
  set(gca,...
    'XLim',sort([w1(1),w1(end)]),...
    'YLim',sort([w3(1),w3(end)]),...
    'ZLim',sort([w5(1),w5(end)]));
  if strcmpi(tick_mode,'man'),
    set(gca,...
      'XTick',ticks,...
      'YTick',ticks,...
      'ZTick',ticks)
  end
  if strcmpi(label_mode,'man'),
    set(gca,...
      'XTickLabel',labels,...
      'YTickLabel',labels,...
      'ZTickLabel',labels);
  end
  daspect([1 1 1]);
  %  view([-45 30]);
  % view([-45 54.142]);
  view(3);
  camproj(projection_type);
  
  %only light the isosurfs
  if n_shells>2
    set([p10],'FaceLighting','gouraud',...
      'DiffuseStrength',0.7,...
      'SpecularStrength',1);
  end
  set([p25],'FaceLighting','gouraud',...
    'DiffuseStrength',0.7,...
    'SpecularStrength',1);
  set([p50],'FaceLighting','gouraud',...
    'DiffuseStrength',0.7,...
    'SpecularStrength',1);
  hl(1) = camlight;
  hl(2) = camlight(-120,45);
  l=line([w1(1),w1(end)],[w1(1),w1(end)],[w1(1),w1(end)]);
  set(l,'Color',[0 0 0]);
  if flag_xlabel,xlabel('\omega_1 / 2\pic');end
  if flag_ylabel,ylabel('\omega_3 / 2\pic');end
  if flag_zlabel,zlabel('\omega_5 / 2\pic');end

  set(gca,'TickDir','out')
  axis vis3d
  
  if nargout>=1
    varargout{1}=ha;
  end

