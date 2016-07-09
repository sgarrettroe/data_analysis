function viewCLSslice(peakFitStruct,ii,jj)

X = peakFitStruct(ii,jj).w3;
Y = peakFitStruct(ii,jj).R;
fitresult = peakFitStruct(ii,jj).fitresult;
XX = X(1):0.1:X(end);
R2 = peakFitStruct(ii,jj).gof.rsquare;
figure,plot(X,Y,'.',XX,fitresult(XX))
title(['R^2 = ',sprintf('%0.5g',R2)])
xlabel('Frequency (cm^{-1})')
ylabel('Intensity (arb. units)')