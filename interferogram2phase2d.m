function [out,q] = interferogram2phase2d(s1,s2,ind)
%convert interferograms to a phase difference
%function out = interferogram2phase2d(s1,s2,pixel);
y1 = s1.signal(:,ind);
y2 = s2.signal(:,ind);
x = s1.time;

dx = x(2)-x(1);
if dx<2.11143
  warning('Time step is a little wrong! Changing to 2.11143 fs/fringe');
  x = x./dx*2.11143;
end

w=fftFreqAxis(x,...
  'fftshift','off',...
  'freq_units','wavenumbers',...
  'time_units','fs',...
  'undersampling',0);

y1hat = fft(ifftshift(y1));
y2hat = fft(ifftshift(y2));
spec1 = abs(y1hat);
spec2 = abs(y2hat);
ph1 = angle(y1hat);
ph2 = angle(y2hat);
dp = ph2-ph1;

[m,i] = max(spec1);
ind2 = (i-4:i+4);
ind3 = (i-1:i+1);
out = mean(dp(ind3));

figure(1),plot(x,y1,x,y2)
figure(2),
h(1)=subplot(2,1,1);
plot(w(ind2),spec1(ind2)./m);
ylabel('spectrum')

h(2)=subplot(2,1,2);
%plot(w(ind2),ph1(ind2),'b-',w(ind2),ph2(ind2),'r-',...
%  w(ind3),ph1(ind3),'bo',w(ind3),ph2(ind3),'ro');
plot(w(ind2),ph1(ind2),'b-',w(ind2),ph2(ind2),'r-')
set(gca,'YLim',[-pi pi]);
xlabel('\omega/2\pic')
ylabel('phase')
linkaxes(h,'x');

q.w = w;
q.spec = spec1;
q.p1 = ph1;
q.p2 = ph2;
q.dp = out;
