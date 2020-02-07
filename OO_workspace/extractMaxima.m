function [out] = extractMaxima(dataobj,options)
% fitStruct = extractMaxima(dataobj,options)
%
% extractMaxima attempts to fit the peak shape of slices of a 2D-IR
% spectrum for use in centerline slope analysis. Input (dataobj) is assumed
% to be (cropped) 2D-IR data, which is further cropped for fitting in this
% function. For each spectrum, this code attempts to fit slices along
% omega_3 with the inputted fitting function (normally a sum of Gaussians
% or Voigt profiles). The function returns the data, fit result, goodness
% of fit, and fitinfo in an n x m structure, where n is the number of t2
% points, and m is the number of omega_1 values fitted.
%
%   required fields for 'options'
%       'range1' - the omega_1 fitting range
%       'range3' - the omega_3 fitting range
%       'fitfcn' - the fitting function
%       'startpoint' - the starting point of the fit 
%       'lb' - the lower bound of the fit
%       'ub' - the upper bound of the fit
%
%   possible additional fields
%       'estimateArea' - (default false) if true, attempts to estimate the
%       area of the peak by calling 'estimatePeakArea.m', and replaces the
%       amplitude (assumed format 'a\d') values in the fitting bounds based
%       on the area. Useful when using a normalized fitting function (like
%       our Voigt profile) 
%
%       'fitPosPeak' - (default false) specifies the behavior of
%       estimatePeakArea. If true, it will attempt to fit a positive peak
%       (probably the ESA). If false, it will fit a negative peak (GSB/SE).
%
%       'flag_plot' - (default false) If true, plots fitting results as
%       they come in.
%
%       'estimate_flag_plot' - (default false) similar to flag_plot, but
%       for estimatePeakArea
%
% TODO: Allow more control over relative magnitudes in estimateArea

% default values
flag_plot = 0;
estimateArea = 0;
fitPosPeak = 0;
estimate_flag_plot = 0;

% checking the fields that the user passed to the function to ensure that
% they match expected fields
required_fields = {'range1' 'range3' 'fitfcn' 'flag_plot' 'lb' 'startpoint' 'ub'};
possible_fields = {'estimateArea' 'estimate_flag_plot' 'fitPosPeak' 'fitfcn' 'flag_plot' 'lb' 'startpoint' 'ub'};
given_fields = fields(options);
missing_fields = setdiff(sort(required_fields),sort(given_fields));
unexpected_fields = setdiff(sort(given_fields),sort(possible_fields));
if ~isempty(missing_fields)
    error('Required options fields not found.');
end
if ~isempty(unexpected_fields)
    warning('Unexpected options fields.')
end

% unpacking the options structrue into a series of variables 
for ji=1:length(given_fields)
    eval([given_fields{ji} '=options.' given_fields{ji} ]);
end

out = struct('fitresult',[],'gof',[],'fitinfo',[]);
for ii = 1:length(dataobj)
    w1 = dataobj(ii).w1;
    w3 = dataobj(ii).w3;
    w3a = min(w3):0.1:max(w3);
    R = dataobj(ii).R;
    for ij = 1:length(w1)
        % redefining our starting parameters for each spectrum slice.
        if estimateArea
            area = estimatePeakArea(w3(:),R(:,ij),'fitPosPeak',fitPosPeak,'flag_plot',estimate_flag_plot);
            % got rid of hard coding -- now amplitudes can be located
            % anywhere in the bounds, and there can be any number of
            % amplitudes.
            coefficients = coeffnames(fitfcn);
            aTest = regexp(coefficients,'a\d'); % find the pattern 'a#'
            ind = ~cellfun(@isempty,aTest); % find where those are living
            startpoint(ind) = area; % redefine our startpoints
            lb(ind) = 0.5*area;
            ub(ind) = 2*area;
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
            if ii == 1 && ij == 1
                fig = figure;
                fig.Color = 'w';
                clf
                p1 = plot(w3(:),R(:,ij),'b.',w3a,fitresult(w3a));
                fig.Children.Box = 'off';
                fig.Children.TickDir = 'out';
                title(sprintf('(%i,%i)',ii,ij));
                xlabel('Frequency (cm^{-1})')
                ylabel('\Delta A (a.u.)')
            else
                p1(1).YData = R(:,ij);
                p1(2).YData = fitresult(w3a);
                fig.Children.Title.String = sprintf('(%i,%i)',ii,ij);
            end
        end
    end
end