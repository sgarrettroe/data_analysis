function [out] = fitEllip(data,options)
startpoint = options.startpoint;
ub = options.upper_bound;
lb = options.lower_bound;
flag_plot = options.flag_plot;

if isfield(options,'fitfcn')
    fitfcn = options.fitfcn;
else
    fitfcn = fittype( @(a, mx, my, sD, sA, delta, c, x, y)...
    a.*(-exp(-(bsxfun(@plus, x-mx, y-my)).^2/(2*sD^2)) .* ...
    exp(-(bsxfun(@minus, x-mx, y-my)).^2/(2*sA^2))...
    + exp(-(bsxfun(@plus, x-mx, y-my+delta)).^2/(2*sD^2)) .* ...
    exp(-(bsxfun(@minus, x-mx, y-my+delta)).^2/(2*sA^2)) ) + c ,...
    'coeff', {'a', 'mx', 'my', 'sD', 'sA', 'delta', 'c'}, ...
    'indep', {'x', 'y'}, 'dep', 'z');
end

out = struct('ellipFit',[]);
for ii = 1:length(data)
    xx = data(ii).w1;
    yy = data(ii).w3;
    zz = data(ii).R;
    [XX,YY] = meshgrid(xx,yy);
    xxx = XX(:);
    yyy = YY(:);
    zzz = zz(:);
    
    A = [-max(max(xx))];
    B = horzcat(A,startpoint);
    [fitresult,gof,info] = fit([xxx,yyy],-zzz,fitfcn,... %%%%
      'StartPoint', B, ...
      'upper',ub ,'lower',lb);

  if flag_plot
    figure(10),clf,my2dPlot(xx,yy,zz,'pumpprobe',false);
    figure(11),clf,my2dPlot(xx,yy,-fitresult(XX,YY),'pumpprobe',false)
  end
  
    out(ii).fitresult = fitresult;
    out(ii).gof = gof;
    out(ii).info = info;
end