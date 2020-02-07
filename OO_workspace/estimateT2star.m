function [T2star,sigmaT2star] = estimateT2star(T2,sigmaT2,T1,sigmaT1,Tor,sigmaTor);

invT2star = 1/T2 - 1/(2*T1) - 1/(3*Tor);
T2star = 1/invT2star;
sigmaInvT2star = sqrt(sigmaT2^2/T2^4 + sigmaT1^2/(4*T1^4) + sigmaTor^2/(9*Tor^4));
sigmaT2star = T2star^2*sigmaInvT2star;