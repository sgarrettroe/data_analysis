function S=petersfft3(V,n)

len = size(V,1);

method = 'linear';
%method = 'spline';
%rows
%disp('rows')
for i = 1:len,
  for j = 1:len
      v = V(j,:,i);
      vi = interp1([v v(end)],linspace(1,len+(n-1)/n,len*n),method);
      %vi = interp1(v,linspace(1,len,len*n));
      w = fft(vi);
      V(j,:,i) = [w(1:len/2) w(end-len/2+1:end)];
  end
end

%columns
%disp('cols')
for i = 1:len,
  for j = 1:len
      v = V(:,j,i);
      vi = interp1([v; v(end)],linspace(1,len+(n-1)/n,len*n),method);
      %vi = interp1(v,linspace(1,len,len*n));
      w = fft(vi);
      V(:,j,i) = [w(1:len/2) w(end-len/2+1:end)];
  end
end

%pages
%disp('pages')
for i = 1:len,
  for j = 1:len
      v = squeeze(V(i,j,:));
      vi = interp1([v; v(end)],linspace(1,len+(n-1)/n,len*n),method);
      %vi = interp1(v,linspace(1,len,len*n));
      w = fft(vi);
      V(i,j,:) = [w(1:len/2) w(end-len/2+1:end)];
  end
end

S=V;
