
%----------------------------------------------------------------------------
%3-point joint-probability distribution (one large figure with
%slices indicated)
%----------------------------------------------------------------------------
%i_fig = i_fig+1;
%figure(i_fig),clf reset
figure(8),clf reset
fontsize=12;
ind = 1;
n_i=length(ind);
count = 0;
for ii = ind
  count = count+1;
  subplot(1,n_i,count)
  v = spec3d{ii};
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
  s = squeeze(v(:,22,:))';
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
  s = squeeze(v(:,11,:))';
  max_s = max(max(s));
  [c,h]=contour(freq,freq,s,[max_s*0.1 max_s*0.25 max_s*0.5]);
  hs2=findobj(h,'type','patch');
  for j = 1:length(hs2)
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
if flag_print
  print(gcf,'-dtiff',[base_name,'_120fs_spec_slices.tif'])
end

%% ----------------------------------------------------------------------------
%slice of joint prob dist +100 cm-1
%----------------------------------------------------------------------------
freq_cut = 100; %freq of the cut

figure(9),clf reset
%ind = [20]; %index of slices
[dummy,ind] = min((freq-freq_cut).^2); %index of slices
ind2=[1:n_t2_array]
n_rows = length(ind2);
n_cols=length(ind);
map = myMapRGB(65);
count = 0;
for i = ind2;
  v = spec3d{i};
  count = count+1;
  for j = 1:n_cols
    subplot(n_rows,n_cols,(count-1)*n_cols+j)
    s = squeeze(v(:,ind(j),:))';
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
    text(-230,180,[num2str(2*t2_array(i)*0.01),' ps'],'Fontsize',10)
  end
end
set(gca,'XTick',[-200 0 200],'XTickLabelMode','auto')
xlabel('\delta\omega_3 / 2\pic','Fontsize',10);
set(gcf,'PaperUnits','inches','PaperPosition',[1 1 2.335 8],'InvertHardcopy','off','Color',[230/255 230/255 1]);
if flag_print
   print(gcf,'-depsc','-adobecset',[base_name,'_p100cm_spec_slices.eps'])
end

%----------------------------------------------------------------------------
%slice of joint prob dist -100 cm-1
%----------------------------------------------------------------------------
figure(10),clf reset
%ind=22;
[dummy,ind] = min((freq+freq_cut).^2); %index of slices
n_rows = length(ind2);
n_cols=length(ind);
map = myMapRGB(65);
count = 0;
for i = ind2;
  v = spec3d{i};
  count = count+1;
  for j = 1:n_cols
    subplot(n_rows,n_cols,(count-1)*n_cols+j)
    s = squeeze(v(:,ind(j),:))';
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
    text(-230,180,[num2str(2*t2_array(i)*0.01),' ps'],'Fontsize',10)
  end
end
set(gca,'XTick',[-200 0 200],'XTickLabelMode','auto')
xlabel('\delta\omega_3 / 2\pic','Fontsize',10);
set(gcf,'PaperUnits','inches','PaperPosition',[1 1 2.335 8],'InvertHardcopy','off','Color',[1 1 230/255]);
if flag_print
   print(gcf,'-depsc','-adobecset',[base_name,'_m100cm_spec_slices.eps'])
end

%% projections of 3d
v = spec3d{1};
s_13 = sum(v,3);
s_35 = squeeze(sum(v,2))';
s_15 = squeeze(sum(v,1))';

figure(1),clf
[a1,a2,a3]=my2dPlot(freq,freq,s_13,11);
axes(a1),
xlabel('\delta\omega_{1} / 2\pic')
ylabel('\delta\omega_{3} / 2\pic')
axes(a2)
line(freq,spec1d./max(spec1d).*max(sum(s_13)),'LineStyle',':')
axes(a3)
line(spec1d./max(spec1d).*max(sum(s_13)),freq,'LineStyle',':')

figure(2),clf
[a1,a2,a3]=my2dPlot(freq,freq,s_35,11);
axes(a1),
xlabel('\delta\omega_{3} / 2\pic')
ylabel('\delta\omega_{5} / 2\pic')
axes(a2)
line(freq,spec1d./max(spec1d).*max(sum(s_35)),'LineStyle',':')
axes(a3)
line(spec1d./max(spec1d).*max(sum(s_35)),freq,'LineStyle',':')

figure(3),clf
[a1,a2,a3]=my2dPlot(freq,freq,s_15,11);
axes(a1),
xlabel('\delta\omega_{1} / 2\pic')
ylabel('\delta\omega_{5} / 2\pic')
axes(a2)
line(freq,spec1d./max(spec1d).*max(sum(s_15)),'LineStyle',':')
axes(a3)
line(spec1d./max(spec1d).*max(sum(s_15)),freq,'LineStyle',':')
