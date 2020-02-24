function [area,min_ind,ymin,FWHM] = estimatePeakArea(x,y,varargin)
% function [area,max_ind,ymax,FWHM] = estimatePeakArea(x,y,varargin)
% [area,max_ind,ymax,FWHM] = estimatePeakArea(x,y) will generate a peak
% area that can be used as input to a Voigt profile (or other normalized
% peak fitting function) to give the the peak fitting code a decent chance
% of having a good starting point and bounds at a variety of w1 values and
% t2 values.
%
% The code was written originally for fitting 2D-IR spectra, and thus it
% assumes a normal 2D-IR peak shape (negative high frequency peak, positive
% low frequency peak). By default it fits the negative peak, though you can
% change this by using the variable argument 'fitPosPeak,' and setting this
% to 'true' (1). Setting this value to false ('0') will fit the negative
% peak, and is the default behavior of the script.
%
% estimatePeakArea(x,y,'fitPosPeak',1)
%
% if you want to see the area you're fitting, use the 'flag_plot' option,
% and set it to '1'.


% default value -- fit the negative peak
fitPosPeak = 0;
flag_plot = 0;
while length(varargin)>=2 %using a named pair
    arg = varargin{1};
    val = varargin{2};
    switch lower(arg)
        case 'fitpospeak'
            fitPosPeak = val;
        case 'flag_plot'
            flag_plot = val;
        otherwise
            warning(['estimatePeakArea: unknown option ',arg])
    end
    varargin = varargin(3:end);
end

[ymin,min_ind] = min(y);

% we're going to interpolate the data to get some extra resolution for
% finding the FWHM.
xx = min(x):0.25:max(x);
yy = spline(x,y,xx);

% A couple of lines to allow us to fit the positive peak instead of the
% negative peak. Although it's counterintuitive, I chose to do this by
% changing the sign of the spectrum and flipping it. This avoids having to
% rewrite a significant chunk of the code.
if fitPosPeak  == true;
    yy = -yy;
    yy = flip(yy);
end

% determine the FWHM of the peak
indices = find(yy <= 0.5.*ymin);
FWHM = xx(max(indices)) - xx(min(indices));


% find the peak minimum and then the region of the peak (on the HF side) in
% the FWHM
[~,yy_ind] = min(yy);
PeakAreaindex = yy_ind:1:yy_ind + length(indices);

% make sure we don't go out of bounds
if max(PeakAreaindex) > numel(xx)
    PeakAreaindex = PeakAreaindex(1):numel(xx);
end
area = abs(2.*trapz(xx(PeakAreaindex),yy(PeakAreaindex)));

if flag_plot
    if fitPosPeak == true;
        a = 1:length(xx);
        a = flip(a);
        ind = find(a >= PeakAreaindex(1) & a <= PeakAreaindex(end));
        yy = flip(yy);
        yy = -yy;
        figure(145);clf
        plot(x,y,'o'),hold on,plot(xx,yy,'k'),plot(xx(ind),yy(ind),'r')
    else
        figure(145);clf,plot(x,y,'o'),hold on,plot(xx,yy,'k'),plot(xx(PeakAreaindex),yy(PeakAreaindex),'r')
    end
    xlabel('x')
    ylabel('y')
    legend('data','interpolation','selected area')
    legend('boxoff')
end
