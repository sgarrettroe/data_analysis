i_fig = 3;
figure(i_fig),clf reset,
if get(gcf,'WindowStyle')=='normal',set(gcf,'Position',[0 394 712 553]);end
iso_value = [0.1 0.25 0.5];
alpha_value = [0.1 0.5 0.75];
%for i=2
for i = 1:n_t2_array
  figure(i_fig);clf
  v = squeeze(spec3d{i}(:,:,:));
  my3dPlot3(freq,freq,freq,v);
  if n_t2_array>1,pause;end
end
orient landscape
if flag_movie,
  fname = [base_name,'_spec3d.avi'];
  make_spec3d_movie
end
