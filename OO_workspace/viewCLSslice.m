function viewCLSslice(peakFitStruct,ii,jj)
% viewCLSslice(peakFitStruct,i,j) Displays the CLS fitting results from
% EXTRACTMAXIMA.
%
% EXTRACTMAXIMA generates a data structure containing the fitresult,
% goodness of fit, fitinfo, and w3, R, and w1 values used for each fitting.
% This code allows you to easily visualize the results.
% 
% The first index 'i' generally corresponds to the indexing within your
% 2D-IR data (i.e. 'i = 1' will corresponds to data(1)). The second index
% 'j' then corresponds to each w1 value in the 2D-IR spectrum for which the
% peaks were fitted.

X = peakFitStruct(ii,jj).w3;
Y = peakFitStruct(ii,jj).R;
fitresult = peakFitStruct(ii,jj).fitresult;
XX = X(1):0.1:X(end);
R2 = peakFitStruct(ii,jj).gof.rsquare;
figure,plot(X,Y,'.',XX,fitresult(XX))
title(['R^2 = ',sprintf('%0.5g',R2)])
xlabel('Frequency (cm^{-1})')
ylabel('Intensity (arb. units)')