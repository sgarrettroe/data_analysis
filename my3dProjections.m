function my3dProjections(freq,spec3d,ind)
%function my3dProjections(freq,spec3d,ind)

ifig = gcf;

% projections of 3d
v = spec3d{ind};
s_13 = sum(v,3);
s_35 = squeeze(sum(v,2))';
s_15 = squeeze(sum(v,1))';

figure(ifig),clf
[a1,a2,a3]=my2dPlot(freq,freq,s_13,10);
axes(a1),
xlabel('\delta\omega_{1} / 2\pic')
ylabel('\delta\omega_{3} / 2\pic')
%axes(a2)
%line(freq,spec1d./max(spec1d).*max(sum(s_13)),'LineStyle',':')
%axes(a3)
%line(spec1d./max(spec1d).*max(sum(s_13)),freq,'LineStyle',':')

figure(ifig+1),clf
[a1,a2,a3]=my2dPlot(freq,freq,s_35,10);
axes(a1),
xlabel('\delta\omega_{3} / 2\pic')
ylabel('\delta\omega_{5} / 2\pic')
%axes(a2)
%line(freq,spec1d./max(spec1d).*max(sum(s_35)),'LineStyle',':')
%axes(a3)
%line(spec1d./max(spec1d).*max(sum(s_35)),freq,'LineStyle',':')

figure(ifig+2),clf
[a1,a2,a3]=my2dPlot(freq,freq,s_15,10);
axes(a1),
xlabel('\delta\omega_{1} / 2\pic')
ylabel('\delta\omega_{5} / 2\pic')
%axes(a2)
%line(freq,spec1d./max(spec1d).*max(sum(s_15)),'LineStyle',':')
%axes(a3)
%line(spec1d./max(spec1d).*max(sum(s_15)),freq,'LineStyle',':')
