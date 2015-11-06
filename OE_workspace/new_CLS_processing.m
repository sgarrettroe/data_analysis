function [peakFit,CLS,CFFit] = new_CLS_processing(data,options,CFoptions)
% NEW_CLS_PROCESSING is for troubleshooting the fitting of our CO2
% data with the centerline slope (CLS) method.
%
% It is assumed that both 'options' and 'CFoptions' will be structures with
% the correct options bundled into them. The general form would be:
%        options.[X] where [X] is the name of the option
% necessary peak fitting options:
%       -range1 (range in w1)
%       -range3 (range in w3)
%       -starting point ('startpoint')
%       -upper bound ('ub')
%       -lower bound ('lb')
%       -lineshape functional form, ('fitfcn') e.g.
%          fitfcn = @(w,center,w_g,w_l,anh,a1,a2) ...
%                    a1.*voigt(w,center,w_g,w_l) - ...
%                    a2.*voigt(w,center - anh,w_g,w_l);
%        For now, you must name your center frequency 'center'. I suggest
%        putting relatively tight bounds on the anharmonicity.
% optional peak fitting options:
%       -flag_plot (a non-zero value will cause the function to plot
%       results as we go)
%       -zero-padding factor (zp_factor) - setting this will change the
%       zero padding on our data, and then re-interpolate it along w1
%
% necessary correlation function options:
%       -starting point ('startpoint')
%       -upper bound ('ub')
%       -lower bound ('lb')
%       -correlation function form ('fitfcn'), e.g.
%           F2Exp = fittype( @(a, a2, t1, t2, c, x) a.*exp(-x./t1)+ ...
%                   a2.*exp(-x./t2)+a3, 'coeff', {'a1', 'a2', 't1', ...
%                   't2', 'a3'}, 'indep', {'x'});
%
% NOT YET IMPLEMENTED: 
%       -- Test for goodness of peak fitting (ie. a Voigt profile
%       does not capture the wings that we see on the experiment), and
%       subsequent limiting of the fitting area


t2_array = [data.t2]./1000;
range1 = options.range1;
range3 = options.range3;
startpoint = options.startpoint;
lb = options.lb;
ub = options.ub;
fitfcn = options.fitfcn;
% how the hell do I avoid demanding that the fitfcn calls its center
% frequency 'center'? Is there a way to dynamically determine this? Bill
% suggested putting an underscore (_) after the center frequency, then
% searching for that.

if isfield(options,'flag_plot')
    flag_plot = options.flag_plot;
else
    flag_plot = 0;
end
CFoptions.flag_plot = flag_plot;
% flag_plot will later determine if we're going to plot our CLS fits
% against the data that generated them.

%let's add some interpolation
if isfield(options,'zp_factor')
    zp_factor = options.zp_factor;
    data = interp2D(data,zp_factor);
end

dataobj = cropData(data,range1,range3);
% We're taking the amount of 'data' that we're going to be working with
% out, so that later we can pass it to other functions more easily and
% manipulate it without having to worry as much about indexing.

peakFit = extractMaxima(dataobj,startpoint,lb,ub,fitfcn,flag_plot);
% This function should generate a structure of fitting results from fitting
% the maximum of the each slice along w3.

[m,n] = size(peakFit);
maxMatrix = zeros(m,n,2);
for ii = 1:m
    for jj = 1:n
         maxMatrix(ii,jj,1) = peakFit(ii,jj).fitresult.center;
         dummy = confint(peakFit(ii,jj).fitresult);
         err_b = dummy(:,1);
         maxMatrix(ii,jj,2) = (err_b(2)-err_b(1))/2;
    end
end
clear dummy err_b m n
% Should populate a matrix of size(data,w1,2) that will contain
% all of the frequencies of the peak maxima in w3 for each element of data
% and each value of w1, and also their standard deviations.

CLS = fitCLS(dataobj,maxMatrix,flag_plot);
% Fitting the CLS based values based on the maxima and standard deviations
% we just calculated.

CFFit = corrFcnFit(t2_array,CLS,CFoptions);

if flag_plot
end