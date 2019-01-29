function S=petersfft(v,n)

len = length(v);
method = 'linear';
vi = interp1([v 0],linspace(1,len+(n-1)/n,len*n),method);
w = fft(vi);
v = [w(1:len/2) w(end-len/2+1:end)];
S=v;
