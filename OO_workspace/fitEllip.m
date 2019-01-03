function [out] = fitEllip(data,options)
% ellipticityStruct = fitEllip(croppedDataStruct,options)
%
% The inputs to this function should be cropped data (over the range you
% want to fit it) including an w1 array, an w3 array, and an R matrix.
% Required options are an ellipticity fitting startpoint, upper bound, and
% lower bound. The form of the fitting function is a bit constrained right
% now.
%
% The basic fitting function used is:
%         fitfcn = fittype( @(a, mx, my, sD, sA, delta, c, x, y)...
%         a.*(-exp(-(bsxfun(@plus, x-mx, y-my)).^2/(2*sD^2)) .* ...
%         exp(-(bsxfun(@minus, x-mx, y-my)).^2/(2*sA^2))...
%         + exp(-(bsxfun(@plus, x-mx, y-my+delta)).^2/(2*sD^2)) .* ...
%         exp(-(bsxfun(@minus, x-mx, y-my+delta)).^2/(2*sA^2)) ) + c ,...
%         'coeff', {'a', 'mx', 'my', 'sD', 'sA', 'delta', 'c'}, ...
%         'indep', {'x', 'y'}, 'dep', 'z');
%
% You can define your own fitting function using an 'options.fitfcn'. With
% how the script is currently written, you want your amplitude to be the
% first element, and sD and sA (diagonal and antidiagonal widths) to be the
% fourth and fifth elements, so that it plays nicely with our later
% scripts.
%
% TODO: Add variable-length input argument list to find the starting
% amplitude with a min or a max, or to not find the starting amplitude.

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

out = struct('fitresult',[],'gof',[],'fitinfo',[]);
for ii = 1:length(data)
    xx = data(ii).w1;
    yy = data(ii).w3;
    zz = data(ii).R;
    [XX,YY] = meshgrid(xx,yy);
    xxx = XX(:);
    yyy = YY(:);
    zzz = zz(:);
    
    A = [-max(max(zz))];
    B = horzcat(A,startpoint);
    [fitresult,gof,fitinfo] = fit([xxx,yyy],-zzz,fitfcn,... %%%%
      'StartPoint', B, ...
      'upper',ub ,'lower',lb);

  if flag_plot
    figure(10),clf,my2dPlot(xx,yy,zz,'pumpprobe',false);
    figure(11),clf,my2dPlot(xx,yy,-fitresult(XX,YY),'pumpprobe',false)
  end
  
    out(ii).fitresult = fitresult;
    out(ii).gof = gof;
    out(ii).fitinfo = fitinfo;
end