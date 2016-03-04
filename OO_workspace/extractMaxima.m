function [out] = extractMaxima(dataobj,options)
% function [out] = extractMaxima(dataobj,startpoint,lb,ub,fitfcn,flag_plot)
startpoint = options.startpoint;
lb = options.lb;
ub = options.ub;
fitfcn = options.fitfcn;
flag_plot = options.flag_plot;

if isfield(options,'estimateArea')
    estimateArea = options.estimateArea;
else
    estimateArea = 0;
end

out = struct('fitresult',[],'gof',[],'fitinfo',[]);
for ii = 1:length(dataobj)
w1 = dataobj(ii).w1;
w3 = dataobj(ii).w3;
w3a = min(w3):0.1:max(w3);
R = dataobj(ii).R;
   for ij = 1:length(w1)
       % redefining our starting parameters for each spectrum slice. I'm
       % pretty sure the Voigt profile I'm using is normalized (ie the
       % total area is 1. So if we can get an estimate of the total area
       % under one of our curves by a trapezoid rule, that would be spiffy.
       if estimateArea
           area = estimatePeakArea(w3(:),R(:,ij));
           startpoint(2:3) = area;
           lb(2:3) = 0.5*area;
           ub(2:3) = 2*area;
       end
       % the fitting part
       [fitresult,gof,fitinfo] = fit(w3(:),R(:,ij),fitfcn, 'StartPoint',startpoint,...
           'lower',lb,'upper',ub);
       out(ii,ij).fitresult = fitresult;
       out(ii,ij).gof = gof;
       out(ii,ij).fitinfo = fitinfo;
       out(ii,ij).w3 = w3;
       out(ii,ij).R = R(:,ij);
       out(ii,ij).w1 = w1(ij);
       if flag_plot
           figure(1),clf,plot(w3(:),R(:,ij),'b.',w3a,fitresult(w3a));
%            pause
       end
     %  out(ii,ij).w1 = w1(ij);
   end
end