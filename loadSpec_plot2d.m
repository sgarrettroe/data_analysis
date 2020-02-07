figure(2),clf
for i=1:n_t2_array,
  my2dPlot4(freq,freq,spec2d{i})
end

% map = myMapRGB(65);
% for i=1:n_t2_array,
%   %subplot(1,n_t2_array,i);
%   figure(i_fig),
%   %contourf(freq,freq,spec2d{i},10);
%   [a1,a2,a3]=my2dPlot(freq,freq,spec2d{i},10);
%   axes(a2)
%   line(freq,spec1d./max(spec1d).*max(sum(spec2d{i})),'LineStyle',':')
%   axes(a3)
%   line(spec1d./max(spec1d).*max(sum(spec2d{i})),freq,'LineStyle',':')
%   %line([freq(1),freq(end)],[freq(1),freq(end)],'Color',[0 0 0])
%   %axis equal, axis tight
%   %colormap(map);
%   %myCaxis(spec2d{i});
%   if(n_t2_array>1),pause,end
% end
