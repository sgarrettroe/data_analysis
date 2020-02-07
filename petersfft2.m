function S=petersfft2(V,n)

len = size(V,1);

method = 'linear';
%rows
%disp('rows')
for j = 1:len
  v = V(j,:);
  vi = interp1([v 0],linspace(1,len+(n-1)/n,len*n),method);
  %  vi = interp1(v,linspace(1,len,len*n));
  w = fft(vi);
  V(j,:) = [w(1:len/2) w(end-len/2+1:end)];
end

%columns
%disp('cols')
for j = 1:len
  v = V(:,j);
  vi = interp1([v; 0],linspace(1,len+(n-1)/n,len*n),method);
  %vi = interp1(v,linspace(1,len,len*n));
  w = fft(vi);
  V(:,j) = [w(1:len/2) w(end-len/2+1:end)];
end

S=V;
