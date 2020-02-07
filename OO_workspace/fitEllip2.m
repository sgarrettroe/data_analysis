function [out] = fitEllip2(dataStruct,coeffStruct,varargin)
% ellipticityStruct = fitEllip2(croppedDataStruct,coeffStruct,varargin)
%
% The inputs to this function should be cropped data (over the range you
% want to fit it) including an w1 array, an w3 array, and an R matrix.
%
%

out = struct('fitresult',[],'gof',[],'fitinfo',[]);

estimateAmplitude = false;
fitfcn = fittype( @(a, mx, my, sD, sA, delta, c, x, y)...
    a.*(-exp(-(bsxfun(@plus, x-mx, y-my)).^2/(2*sD^2)) .* ...
    exp(-(bsxfun(@minus, x-mx, y-my)).^2/(2*sA^2))...
    + exp(-(bsxfun(@plus, x-mx, y-my+delta)).^2/(2*sD^2)) .* ...
    exp(-(bsxfun(@minus, x-mx, y-my+delta)).^2/(2*sA^2)) ) + c ,...
    'coeff', {'a', 'mx', 'my', 'sD', 'sA', 'delta', 'c'}, ...
    'indep', {'x', 'y'}, 'dep', 'z');
flag_plot = 0;
ampString = 'a';

startPoint = coeffStruct.startpoint;
lb = coeffStruct.lb;
ub = coeffStruct.ub;

% Deal with variable user inputs
while length(varargin)>=2 %using a named pair
    arg = varargin{1};
    val = varargin{2};
    switch lower(arg)
        case 'estimateamplitude'
            estimateAmplitude = true;
            if strcmp(val,'max')||strcmp(val,'min')
                ampType = val;
            else
                warning(['Unknown method for estimating amplitude: ',val,...
                    '. Allowed values are: "max" or "min. Not estimating amplitudes.'])
                estimateAmplitude = false;
            end
        case 'flagplot'
            flag_plot = val;
        case 'fitfcn'
            fitfcn = val;
        case 'amplitudestring'
            ampString = val;
    end
    varargin = varargin(3:end);
end

coeffNameArray = coeffnames(fitfcn);

for ii = 1:length(dataStruct)
    x = dataStruct(ii).w1;
    y = dataStruct(ii).w3;
    z = dataStruct(ii).R;    
    
    if estimateAmplitude == true
        ampIndex = strcmp(coeffNameArray,ampString);
        switch lower(ampType)
            case 'max'
                ampEstimate = -max(max(z));
            case 'min'
                ampEstimate = min(min(z));
        end
        startPoint(ampIndex) = ampEstimate;
    end

    [xx,yy,zz] = prepareSurfaceData(x,y,z);  
    [fitresult,gof,fitinfo] = fit([xx,yy],-zz,fitfcn,... %%%%
      'StartPoint', startPoint, ...
      'upper',ub ,'lower',lb);

  if flag_plot
    figure(10),clf,my2dPlot(x,y,z,'pumpprobe',false);
    [X,Y] = meshgrid(x,y);
    figure(11),clf,my2dPlot(x,y,-fitresult(X,Y),'pumpprobe',false)
  end
  
    out(ii).fitresult = fitresult;
    out(ii).gof = gof;
    out(ii).fitinfo = fitinfo;
end