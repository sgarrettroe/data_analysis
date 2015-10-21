%% cd to the wherever the data is and load the data
clear data out out_wt ellip_out c2
process2dDataInitialize
cd(localdir{1});

% load designated _summary.mat file
file_name = sprintf('%s-%02i_summary.mat',datestring{1},experiment{1});
load(file_name);

process2dDataZhesFcns

ellipticity.ind1 = find(data(1).w1>ellipticity.range1(1) & data(1).w1<ellipticity.range1(2));
ellipticity.ind3 = find(data(1).w3>ellipticity.range3(1) & data(1).w3<ellipticity.range3(2));


%% #11 Generate ellipticities
% Initialize variables used in the loop
std_sD = zeros(size(data));
std_sA = zeros(size(data));
sd = zeros(size(data));
sa = zeros(size(data));
anharmonicity = zeros(size(data));
de_sD = zeros(size(data)); %derivatives for sD
de_sA = zeros(size(data)); %derivatives for sA
std_e = zeros(size(data));

% Fit for ellipticity
for ii=1:length(data); %thisscanlist{1}
  xx = data(ii).w1(ellipticity.ind1);
  yy = data(ii).w3(ellipticity.ind3); 
  zz = data(ii).R(ellipticity.ind3,ellipticity.ind1);
  [XX,YY] = meshgrid(xx,yy);
  xxx = XX(:);
  yyy = YY(:);
  zzz = zz(:);
  
%     F2DS = fittype( @(a, mx, my, sD, sA, delta, c, x, y)...
%       a.*(-exp(-(bsxfun(@plus, x-mx, y-my)).^2/(2*sD^2)) .* ...
%       exp(-(bsxfun(@minus, x-mx, y-my)).^2/(2*sA^2))...
%       + exp(-(bsxfun(@plus, x-mx, y-my+delta)).^2/(2*sD^2)) ...
%       .* exp(-(bsxfun(@minus, x-mx, y-my+delta)).^2/(2*sA^2)) )...
%       +c,...
%       'coeff', {'a', 'mx', 'my', 'sD', 'sA', 'delta', 'c'}, ...
%       'indep', {'x', 'y'}, 'dep', 'z');

    A = [-max(max(xx))];
    B = horzcat(A,ellipticity.startpoint);
    out = fit([xxx,yyy],-zzz,F2DS,...
      'StartPoint', B, ...
      'upper',ellipticity.upper_bound ,... 
      'lower',ellipticity.lower_bound);

  figure(10),clf,my2dPlot(xx,yy,zz,'pumpprobe',false);
  figure(11),clf,my2dPlot(xx,yy,-out(XX,YY),'pumpprobe',false)
  
  data(ii).ellipticity.fit = out;
  %get amplitude here out.a?
  dummy = confint(data(ii).ellipticity.fit);
  err_sD = dummy(:,4);
  err_sA = dummy(:,5);
  std_sD(ii) = (err_sD(2) - err_sD(1))/2; %std of sd
  std_sA(ii) = (err_sA(2) - err_sA(1))/2; %std of sa
  sd(ii) = out.sD;
  sa(ii) = out.sA;
  sD = out.sD;
  sA = out.sA;
  anharmonicity(ii) = out.delta;
  de_sD(ii) = (4*sD*sA^2)/(sD^2+sA^2)^2; %derivatives for sD
  de_sA(ii) = -(4*sD^2*sA)/(sD^2+sA^2)^2; %derivatives for sA
  std_e(ii) = (de_sD(ii).^2*std_sD(ii).^2 + de_sA(ii).^2*std_sA(ii).^2)^0.5;
end
 ellipticity.ellip=(sd.^2-sa.^2)./(sd.^2+sa.^2);
 ellipticity.t = [data.t2];
 ellipticity.std_e = std_e;

% plot data with errorbar
% figure(99),errorbar(ellipticity.t,ellipticity.ellip,ellipticity.std_e,'rx');

%% Fitting the Ellipticity to a Functional Form
% Assign and normalize weights
% weight is assigned to be 1 over std squared
clearvars wt
wt = 1./(ellipticity.std_e.^2);
wt_mean = mean(wt);
wt = wt/wt_mean;

% fitting type selection
% F2Exp = fittype( @(a, b, t1, t2, c, x) a.*exp(-x./t1)+b.*exp(-x./t2)+c, 'coeff', {'a', 'b', 't1', 't2', 'c'}, 'indep', {'x'});
% F2Exp_nc = fittype( @(a, b, t1, t2, x) a.*exp(-x./t1)+b.*exp(-x./t2), 'coeff', {'a', 'b', 't1', 't2'}, 'indep', {'x'});

% fit with weight included
% origin fit
clearvars out out_wt

% weighted fit
ellipticity.fit.biexp_const = fit(ellipticity.t(:), ellipticity.ellip(:), ...
    F2Exp, 'StartPoint', [1, 1, 2e3, 2e4,0],'lower', [0, 0, 0, 0,0],'Weight',wt);
ellipticity.fit.biexp = fit(ellipticity.t(:), ellipticity.ellip(:), F2Exp_nc, ...
    'StartPoint', [1, 1, 2e3, 2e4],'lower', [0, 0, 0, 0],'Weight',wt);
%%
A = max(ellipticity.t);
t = ellipticity.t(1):A/100:A;

figure(98),clf
    hold on    
    errorbar(ellipticity.t,ellipticity.ellip,ellipticity.std_e,'rx')
    plot(t, ellipticity.fit.biexp_const(t),'r--')
    plot(t, ellipticity.fit.biexp(t),'b--')
    hold off
    h = gca;
    set(h,'Tickdir','out')
    legend('Data','Biexponential + Const','Biexponential')
    ylabel('Ellipticity (c_2)')
    xlabel('t_2 / fs')
    title(data(1).name)
clear A t h
%%
beep
flag_save = input('Save output? (1 = "YES"; 0 = "NO"): ');
if flag_save
    save(sprintf('%s-%02i_summary.mat',datestring{dirindex(1)},experiment{1}),...
        'data','ellipticity','-append')
end
flag_save = 0;