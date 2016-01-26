function [area,min_ind,ymin,FWHM] = estimatePeakArea(x,y)
% function [area,max_ind,ymax,FWHM] = estimatePeakArea(x,y)
% [area,max_ind,ymax,FWHM] = estimatePeakArea(x,y) will generate a peak
% area that can be used as input to a Voigt profile (or other normalized
% peak fitting function) to give the the peak fitting code a decent chance
% of having a good starting point and bounds at a variety of w1 values and
% t2 values.
%
% This code is written for the CLS fitting of 2D-IR spectra. It assumes
% that you have a simple peak structure (not multiple peaks rolled
% together.
%
% It's currently set up to do a cubic spline interpolation to get some
% extra 'resolution' so I have a better idea of what our FWHM is -- for CO2
% we only have a few points otherwise, and it's too hit or miss if one
% slightly above or below the half maximum.
%

% [ymax,max_ind] = max(y);
[ymin,min_ind] = min(y);

% we're going to interpolate the data to get some extra resolution for
% finding the FWHM. Don't worry, we're not going to use this fit -- just
% use it to define a range to do a trapezoid rule over.
xx = min(x):0.25:max(x);
yy = spline(x,y,xx);

% indices = find(yy >= 0.5.*ymax);
indices = find(yy <= 0.5.*ymin);
FWHM = xx(max(indices)) - xx(min(indices));

% [~,yy_ind] = max(yy);
[~,yy_ind] = min(yy);
PeakAreaindex = yy_ind:1:yy_ind + length(indices);

% area = 2.*trapz(xx(PeakAreaindex),yy(PeakAreaindex));
area = abs(2.*trapz(xx(PeakAreaindex),yy(PeakAreaindex)));