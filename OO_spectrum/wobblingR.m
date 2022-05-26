function [R] = wobblingR(C,stark_effect_order)
% wobbling in a cone model of Kramer Nishda Fayer JCP 2016

switch stark_effect_order
    case 1 % from Kramer2015a
R.para = @(p,x)(3/25)*((11.*C{1}(p,x) + 4*C{3}(p,x))./(1 + 0.8*C{2}(p,x)));
R.perp = @(p,x)(3/25)*(( 7.*C{1}(p,x) - 2*C{3}(p,x))./(1 - 0.4*C{2}(p,x)));
R.iso = @(p,x)C{1}(p,x);
    case 2 % from 
R.para = @(p,x)(1/175)*((28 + 215.*C{2}(p,x) + 72*C{4}(p,x))./(1 + 0.8*C{2}(p,x)));
R.perp = @(p,x)(1/175)*((-14 + 155.*C{2}(p,x) - 36*C{4}(p,x))./(1 - 0.4*C{2}(p,x)));
R.iso = @(p,x)C{2}(p,x);
    otherwise
        error('only order 1 and 2 coded');
end

%orientation
%oPara = @(tau_o,t2) (1/9).*(1+4/5.*exp(-t2./tau_o));
%oPerp = @(tau_o,t2) (1/9).*(1-2/5.*exp(-t2./tau_o));



