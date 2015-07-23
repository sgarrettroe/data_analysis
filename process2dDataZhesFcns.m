%% process2dDataZhesFcns

% Functions for fitting with ellipticity and centerline slope methods

% #01 Generate pure 1D Gausian
SamplePoints = @(mean) [mean-5:0.1:mean+5];
Data1DGclean = @(xx, sigma) exp(-(xx-mean(xx)).^2/(2*sigma^2));
Data1DGnoisy = @(xx, sigma, snr) awgn(Data1DGclean(xx, sigma), snr);
Data2DGclean = @(x1, sigma1, x2, sigma2) exp(-x1.^2/(2*sigma1^2))' * exp(-x2.^2/(2*sigma2^2));
Data2DGnoisy = @(x1, sigma1, x2, sigma2, snr) awgn(Data2DGclean(x1, sigma1, x2, sigma2), snr);
Data1DSclean = @(xx, sigma, delta) exp(-(xx-mean(xx)).^2/(2*sigma^2))-exp(-(xx-mean(xx)-delta).^2/(2*sigma^2));
Data1DSnoisy = @(xx, sigma, delta, snr) awgn(Data1DSclean(xx, sigma, delta), snr);
Data2DSclean = @(xx, sigmaD, yy, sigmaA, delta) -exp(-(bsxfun(@plus, xx', yy)).^2/(2*sigmaD^2)) .* exp(-(bsxfun(@minus, xx', yy)).^2/(2*sigmaA^2)) + exp(-(bsxfun(@plus, xx', yy+delta)).^2/(2*sigmaD^2))' .* exp(-(bsxfun(@minus, xx', yy+delta)).^2/(2*sigmaA^2))';
Data2DSnoisy = @(xx, sigmaD, yy, sigmaA, delta, snr) awgn(Data2DSclean(xx, sigmaD, yy, sigmaA, delta), snr);

ZheFit = @(t, a1, t1, a2, t2, c ) a1^2 * exp(-t/t1) + a2^2 * exp(-t/t2)+c^2;

% #04 Define fit types
F1DT = fittype( @(a,m,s,c,x) a.*exp(-(x-m).^2/(2*s^2))+c, 'coeff',{'a','m','s','c'},'indep',{'x'});%1D
%F2DT = fittype( @(a,s1,s2,c,x,y) a.*exp(-x.^2/(2*s1^2))' * exp(-y.^2/(2*s2^2))+c, 'coeff',{'a','s1','s2','c'},'indep',{'x','y'},'dep','z');%2D
F2DT = fittype( @(a,s1,s2,c,x,y) a.*exp(-x.^2/(2*s1^2)) .* exp(-y.^2/(2*s2^2))+c, 'coeff',{'a','s1','s2','c'},'indep',{'x','y'},'dep','z');%2D
F1DS = fittype( @(a, m, s, delta, c, x) a.*(exp(-(x-m).^2/(2*s^2)) - exp(-(x-m-delta).^2/(2*s^2)))+c, 'coeff',{'a','m','s','delta','c'}, 'indep',{'x'});
%@@@ Add separate mean for 
F2DS = fittype( @(a, mx, my, sD, sA, delta, c, x, y)...
  a.*(-exp(-(bsxfun(@plus, x-mx, y-my)).^2/(2*sD^2)) .* exp(-(bsxfun(@minus, x-mx, y-my)).^2/(2*sA^2))...
  + exp(-(bsxfun(@plus, x-mx, y-my+delta)).^2/(2*sD^2)) .* exp(-(bsxfun(@minus, x-mx, y-my+delta)).^2/(2*sA^2)) )...
  +c,...
  'coeff', {'a', 'mx', 'my', 'sD', 'sA', 'delta', 'c'}, 'indep', {'x', 'y'}, 'dep', 'z');

FExp = fittype( @(a, t, c, x) a.*exp(-x./t)+c, 'coeff', {'a', 't', 'c'}, 'indep', {'x'});
F2Exp = fittype( @(a, b, t1, t2, c, x) a.*exp(-x./t1)+b.*exp(-x./t2)+c, 'coeff', {'a', 'b', 't1', 't2', 'c'}, 'indep', {'x'});
F2Exp_nc = fittype( @(a, b, t1, t2, x) a.*exp(-x./t1)+b.*exp(-x./t2), 'coeff', {'a', 'b', 't1', 't2'}, 'indep', {'x'}); % No constant
F1DP = fittype( @(a, b, c, x) a.*x.^2+b.*x+c, 'coeff', {'a','b','c'}, 'indep', {'x'});
F1DC = fittype( @(a, b, c, d, x) a.*x.^3+b.*x.^2+c.*x+d, 'coeff', {'a','b','c','d'}, 'indep', {'x'});
d = 25;
d1 = 23.8;
FGaus_p = fittype( @(a, b, c, x) -a.*exp(-((x-b)./c).^2), 'coeff', {'a',  'b',  'c'}, 'indep', {'x'});
F2Gaus_p = fittype( @(a, b, c, x) -a.*exp(-((x-b)./c).^2) + ...
    a.*exp(-((x-b+d)./c).^2),...
    'coeff', {'a',  'b',  'c'}, 'indep', {'x'});
F2GausCO2_p = fittype( @(a, b, c, x) -a.*exp(-((x-b)./c).^2) + ...
    a.*exp(-((x-b+d1)./c).^2),...
    'coeff', {'a',  'b',  'c'}, 'indep', {'x'});
F2Lor_p = fittype( @(a, b, c, x)  -a.*c./((x-b).^2 + c.^2) ...
    + a.*c./((x-b+d1).^2 + c.^2), ...
    'coeff', {'a','b','c'},'indep',{'x'}); 