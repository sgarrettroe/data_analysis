function myRPlots(s,ipix)
%plot Rs in time domain at given pixel
n_contours = 10;
map = myMapRGB(n_contours);

subplot(2,3,1)
contourf(s.time,s.time,s.R1(:,:,ipix),n_contours)
axis square
colormap(map)
%myCaxis2(squeeze(s.R1(:,:,ipix)),n_contours)

subplot(2,3,2)
contourf(s.time,s.time,s.R2(:,:,ipix),n_contours)
axis square
colormap(map)
%myCaxis2(squeeze(s.R2(:,:,ipix)),n_contours)
%myCaxis(s.R2(:,:,ipix))

subplot(2,3,4)
contourf(s.time,s.time,s.R3(:,:,ipix),n_contours)
axis square
colormap(map)
%myCaxis2(squeeze(s.R3(:,:,ipix)),n_contours)
%myCaxis(s.R3(:,:,ipix))

subplot(2,3,5)
contourf(s.time,s.time,s.R4(:,:,ipix),n_contours)
axis square
colormap(map)
%myCaxis2(squeeze(s.R4(:,:,ipix)),n_contours)
%myCaxis(s.R4(:,:,ipix))

subplot(2,3,3)
plot(s.time,s.R1(1,:,ipix),s.time,s.R2(1,:,ipix),...
  s.time,s.R3(1,:,ipix),s.time,s.R4(1,:,ipix))

%subplot(2,3,6)
%plot(s.time,s.R1(:,1,ipix),s.time,s.R2(:,1,ipix),...
%  s.time,s.R3(:,1,ipix),s.time,s.R4(:,1,ipix))
subplot(2,3,6)
plot(s.freq,squeeze(s.R1(1,1,:)))
