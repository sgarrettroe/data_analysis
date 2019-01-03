function [norm,s] = correctAmpDrift3d(s)
%correctAmpDrift3d.m try to correct amplitude drift during a measurement by
%normalizing data to the t1=t3=0 fs dispersed data
%call as 
%         [norm,s] = correctAmpDrift3d(s)
% where s is an array of scan structures ,for example a typical use with
% correctPhaseDrift3d would be: 
% 
% dph = zeros(1,9); %phase array
% [norm,test(3:9)] = correctAmpDrift3d(hod(3:9))
% dph(3:9) = correctPhaseDrift3d(hod(3:9));


if ~isstruct(s),
  error('input needs to be a scan structure.')
end

freq = s(1).freq; %spectrometer axis
n_w5 = length(freq);

n_scans = length(s);
uh = cell(1,n_scans);
norm = cell(1,n_scans);
offset = cell(1,n_scans);
ref_scan = 1; 

%loop over diagrams comparing the size (p-p) of the echo signal
%iscan=1;
for iscan = 1:n_scans
  norm{iscan} = zeros(1,4);
  offset{iscan} = zeros(1,4);
  
  %grab the t1=t3=0 spectrum
  uh{iscan}(1,1:n_w5) = squeeze(s(iscan).R1(1,1,:));
  uh{iscan}(2,1:n_w5) = squeeze(s(iscan).R2(1,1,:));
  uh{iscan}(3,1:n_w5) = squeeze(s(iscan).R3(1,1,:));
  uh{iscan}(4,1:31) = squeeze(s(iscan).R4(1,1,:));
  
  
  norm{iscan}(1) = max(uh{iscan}(1,:))-min(uh{iscan}(1,:));
  norm{iscan}(2) = max(uh{iscan}(2,:))-min(uh{iscan}(2,:));
  norm{iscan}(3) = max(uh{iscan}(3,:))-min(uh{iscan}(3,:));
  norm{iscan}(4) = max(uh{iscan}(4,:))-min(uh{iscan}(4,:));
  offset{iscan}(1) = mean(uh{iscan}(1,:));
  offset{iscan}(2) = mean(uh{iscan}(2,:));
  offset{iscan}(3) = mean(uh{iscan}(3,:));
  offset{iscan}(4) = mean(uh{iscan}(4,:));
  if iscan == ref_scan
    ref = norm{iscan}(1);
    ref_offset = offset{iscan}(1);
  end
  norm{iscan} = norm{iscan}./ref;
  offset{iscan} = offset{iscan}-ref_offset;
  for i = 1:4
    uh{iscan}(i,:) = (uh{iscan}(i,:)-offset{iscan}(i))./norm{iscan}(i);
  end
  
  figure(11),clf,plot(freq,uh{iscan}(1,:),freq,uh{iscan}(2,:),freq,uh{iscan}(3,:),freq,uh{iscan}(4,:))
  if n_scans>1,pause,end
end

% try to correct the amplitudes now
for i = 1:n_scans
  s(i).R1 = s(i).R1./norm{i}(1);
  s(i).R2 = s(i).R2./norm{i}(2);
  s(i).R3 = s(i).R3./norm{i}(3);
  s(i).R4 = s(i).R4./norm{i}(4);
end


figure(11),clf
hold on
for i = 1:n_scans
  plot(freq,uh{i}(1,:),...
  freq,uh{i}(2,:),...
  freq,uh{i}(3,:),...
  freq,uh{i}(4,:))
end
hold off

