function dph = correctPhaseDrift3d(s);
%correctPhaseDrift3d.m Try to correct for phase drift between scans by
%looking at the R1 interferogram's phase in the middle of the dispersed
%spectrum (hope that it is enough). Use the dph to correct the inputs to
%absorptive3d.m

if ~isstruct(s) 
  error('input must be a scan structure');
end

%for plotting
t = s(1).time;
w = fftFreqAxis(t,'time_units','fs','fftshift','off','undersampling',s(1).undersampling)

%initialize vars
n_scans = length(s);
uhh = cell(1,n_scans); %time domain
uhh_hat = cell(1,n_scans); %freq domain
mag = cell(1,n_scans);
phase = cell(1,n_scans);
ref_scan = 1; 
ref_pixel = 12;

%grab the time domain data
for iscan = 1:n_scans
  uhh{iscan} = squeeze(s(iscan).R1(1,:,ref_pixel));
end

%do the fft
for i = 1:n_scans
  uhh_hat{i} = fft(uhh{i});
  mag{i} = abs(uhh_hat{i});
  phase{i} = unwrap(angle(uhh_hat{i}));
end

%find the peak of the spectrum
[dummy,imax] = max(mag{ref_scan});%index of peak
disp(['using max at ' num2str(w(imax)) ' cm-1']);

%calculate the phase differences
for i = 1:n_scans
  dph(i) = phase{i}(imax)-phase{1}(imax);
  dph(i) = mod(dph(i),2*pi); %reduce to range [0 360)
  disp(['dph = ' num2str(dph(i)*180/pi) ' deg']);
end

figure(12),
subplot(4,1,1),plot(t,uhh{1},t,uhh{2})
subplot(4,1,2),plot(w,mag{1},w,mag{2})
subplot(4,1,3),plot(w,phase{1},w,phase{2})
subplot(4,1,4),plot(w,real(uhh_hat{1}.*exp(1i*dph(1))),...
  w,real(uhh_hat{2}.*exp(1i*dph(2))))

