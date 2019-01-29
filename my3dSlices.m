function my3dSlices(freq,spec3d,t2_array,t4_array,ind,ind1,ind2)
%function my3dSlices(freq,spec3d,t2_array,t4_array,ind,ind1,ind2)
%
% this takes a 3d spectrum or probability distribution and plots the 3d
% object with slices indicated in one figure, and a time series of the
% slices in another figure
%
% see also my3dSlices2 for just plotting the cuts of a single spectrum or
% probability distribution

map = myMapRGB(65);
i_fig = 100;

%% slices of 3D
%----------------------------------------------------------------------------
%3-point joint-probability distribution (one large figure with
%slices indicated)
%----------------------------------------------------------------------------
i_fig = i_fig+1;
figure(i_fig),clf reset
%figure(8),clf reset
fontsize=12;

n_i=length(ind);
count = 0;
for i = ind
  count = count+1;
  subplot(1,n_i,count)
  v = spec3d{i};
  v=v./max(max(max(v)));
  hold on
  %xy plane
  s=sum(v,3);
  [c,h]=contourf(freq,freq,s,6);
  h=findobj(h,'type','patch');
  colormap(map),myCaxis(s);
  for i=1:length(h),zd=freq(1).*ones(size(get(h(i),'XData')));set(h(i),'ZData',zd);end
  set(h,'FaceLighting','none');
  
  %yz plane
  s=squeeze(sum(v,2))';
  [c,h]=contourf(freq,freq,s,6);colormap(map),myCaxis(s);
  h=findobj(h,'type','patch');
  for i=1:length(h),
    xd=freq(end).*ones(size(get(h(i),'YData')));
    zd=get(h(i),'YData');
    yd=get(h(i),'XData');
    set(h(i),'XData',xd,'YData',yd,'ZData',zd);
  end
  set(h,'FaceLighting','none');
  
  %xz plane
  s=squeeze(sum(v,1))';
  [c,h]=contourf(freq,freq,s,6);colormap(map),myCaxis(s);
  h=findobj(h,'type','patch');
  for i=1:length(h),
    xd=get(h(i),'XData');
    zd=get(h(i),'YData');
    yd=freq(end).*ones(size(get(h(i),'YData')));
    set(h(i),'XData',xd,'YData',yd,'ZData',zd);
  end
  set(h,'FaceLighting','none');
  
  %cut1
  s = squeeze(v(:,ind2,:))';
  max_s = max(max(s));
  [c,h]=contour(freq,freq,s,[max_s*0.1 max_s*0.25 max_s*0.5]);
  hs=findobj(h,'type','patch');
  for j = 1:length(hs)
    xd = freq(22)*ones(size(get(hs(j),'XData')));
    yd = get(hs(j),'XData');
    zd = get(hs(j),'YData');
    set(hs(j),'XData',xd,'YData',yd,'ZData',zd);
  end
  set(hs,'EdgeColor',[0 0 1],'LineWidth',2,'EdgeLighting','none')
  
  %cut2
  s = squeeze(v(:,ind1,:))';
  max_s = max(max(s));
  [c,h]=contour(freq,freq,s,[max_s*0.1 max_s*0.25 max_s*0.5]);
  hs2=findobj(h,'type','patch');
  for j = 1:length(hs)
    xd = freq(11)*ones(size(get(hs2(j),'XData')));
    yd = get(hs2(j),'XData');
    zd = get(hs2(j),'YData');
    set(hs2(j),'XData',xd,'YData',yd,'ZData',zd);
  end
  set(hs2,'EdgeColor',[1 1 0],'LineWidth',2,'EdgeLighting','none')
  
  p10 = patch(isosurface(freq,freq,freq,v,0.10));
  p25 = patch(isosurface(freq,freq,freq,v,0.25));
  p50 = patch(isosurface(freq,freq,freq,v,0.5));
  isonormals(freq,freq,freq,v,p10);
  set(p10,'FaceColor','red','EdgeColor','none','FaceAlpha',0.1);
  set(p25,'FaceColor','red','EdgeColor','none','FaceAlpha',0.5);
  set(p50,'FaceColor','red','EdgeColor','none','FaceAlpha',0.75);
  set(gca,'XLim',[freq(1),freq(end)],...
    'YLim',[freq(1),freq(end)],...
    'ZLim',[freq(1),freq(end)],...
    'XTick',[-200 0 200],...
    'YTick',[-200 0 200],...
    'ZTick',[-200 0 200],...
    'FontSize',fontsize);
  l=line([freq(1),freq(end)],[freq(1),freq(end)],[freq(1),freq(end)]);
  set(l,'Color',[0 0 0]);
  xlabel('\delta\omega_1 / 2\pic','Fontsize',fontsize);ylabel('\delta\omega_3 / 2\pic','Fontsize',fontsize);zlabel('\delta\omega_5 / 2\pic','Fontsize',fontsize);
  daspect([1 1 1]);
  view(3),box off
  %only light the isosurfs
  set([p10 p25 p50],'FaceLighting','gouraud',...
    'DiffuseStrength',0.7,...
    'SpecularStrength',1);
  hl(1) = camlight;
  hl(2) = camlight(-120,45);
 
  set(hs,'EdgeLighting','none')
  set(hs2,'EdgeLighting','none')
  hold off
end
set(gcf,'PaperUnits','inches','PaperPosition',[1 1 6 6]);
%if flag_print
%  print(gcf,'-dtiff',[base_name,'_120fs_spec_slices.tif'])
%end

%% ----------------------------------------------------------------------------
%slice of joint prob dist +100 cm-1
%----------------------------------------------------------------------------
i_fig = i_fig+1;
figure(i_fig),clf reset
%figure(9),clf reset
n_rows = length(ind);
n_cols=length(ind2);
map = myMapRGB(65);
count = 0;
for i = ind;
  v = spec3d{i};
  count = count+1;
  for j = 1:n_cols
    subplot(n_rows,n_cols,(count-1)*n_cols+j)
    s = squeeze(v(:,ind2(j),:))';
    contourf(freq,freq,s,6);
    set(gca,'XLim',[freq(1),freq(end)],...
	    'YLim',[freq(1),freq(end)],...
	    'XTick',[-200 0 200],...
	    'YTick',[-200 0 200],...
	    'XTickLabel',[]);
    l=line([freq(1),freq(end)],[freq(1),freq(end)]);
    set(l,'Color',[0 0 0]);
    daspect([1 1 1]);
    colormap(map);
    myCaxis(s);
    set(gca,'FontSize',8,'Tickdir','out','TickLength',[0.05 0.1]);
    ylabel(['\delta\omega_5 / 2\pic'],'Fontsize',10);
    text(-280,180,[num2str(t2_array(i)+t4_array(i)),' ps'],'Fontsize',10)
  end
end
set(gca,'XTick',[-200 0 200],'XTickLabelMode','auto')
xlabel('\delta\omega_3 / 2\pic','Fontsize',10);
set(gcf,'PaperUnits','inches','PaperPosition',[1 1 2.335 8],'InvertHardcopy','off','Color',[230/255 230/255 1]);
%if flag_print
%   print(gcf,'-depsc','-adobecset',[base_name,'_p100cm_spec_slices.eps'])
%end

%----------------------------------------------------------------------------
%slice of joint prob dist -100 cm-1
%----------------------------------------------------------------------------
i_fig = i_fig+1;
figure(i_fig),clf reset
%figure(10),clf reset
n_rows = length(ind);
n_cols=length(ind1);
map = myMapRGB(65);
count = 0;
for i = ind;
  v = spec3d{i};
  count = count+1;
  for j = 1:n_cols
    subplot(n_rows,n_cols,(count-1)*n_cols+j)
    s = squeeze(v(:,ind1(j),:))';
    contourf(freq,freq,s,6);
    set(gca,'XLim',[freq(1),freq(end)],...
	    'YLim',[freq(1),freq(end)],...
	    'XTick',[-200 0 200],...
	    'YTick',[-200 0 200],...
	    'XTickLabel',[]);
    l=line([freq(1),freq(end)],[freq(1),freq(end)]);
    set(l,'Color',[0 0 0]);
    daspect([1 1 1]);
    colormap(map);
    myCaxis(s);
    set(gca,'FontSize',8,'Tickdir','out','TickLength',[0.05 0.1]);
    ylabel(['\delta\omega_5 / 2\pic'],'Fontsize',10);
    text(-280,180,[num2str(t2_array(i)+t4_array(i)),' ps'],'Fontsize',10)
  end
end
set(gca,'XTick',[-200 0 200],'XTickLabelMode','auto')
xlabel('\delta\omega_3 / 2\pic','Fontsize',10);
set(gcf,'PaperUnits','inches','PaperPosition',[1 1 2.335 8],'InvertHardcopy','off','Color',[1 1 230/255]);
%if flag_print
%   print(gcf,'-depsc','-adobecset',[base_name,'_m100cm_spec_slices.eps'])
%end
