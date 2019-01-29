function s = removeSpectrum3d(s) 
%remove the constant background from the 3d data so I can zeropad without
%that weird artifact from the step function
for i = 1:length(s.freq)
  m = mean(mean(s.R1(:,:,i)));
  s.R1(:,:,i) = s.R1(:,:,i)-m;

  m = mean(mean(s.R2(:,:,i)));
  s.R2(:,:,i) = s.R2(:,:,i)-m;

  m = mean(mean(s.R3(:,:,i)));
  s.R3(:,:,i) = s.R3(:,:,i)-m;

  m = mean(mean(s.R4(:,:,i)));
  s.R4(:,:,i) = s.R4(:,:,i)-m;
end
