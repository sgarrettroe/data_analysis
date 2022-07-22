function [R] = wobblingRv2(C,stark_effect_order)
% wobbling in a cone model of Kramer Nishda Fayer JCP 2016

switch stark_effect_order
    case 1 % from Kramer2015a
R.para = @(x,p)(3/25)*((11.*C{1}(x,p) + 4*C{3}(x,p))./(1 + 0.8*C{2}(x,p)));
R.perp = @(x,p)(3/25)*(( 7.*C{1}(x,p) - 2*C{3}(x,p))./(1 - 0.4*C{2}(x,p)));
R.iso = @(x,p)C{1}(x,p);
    case 2 % from 
R.para = @(x,p)(1/175)*((28 + 215.*C{2}(x,p) + 72*C{4}(x,p))./(1 + 0.8*C{2}(x,p)));
R.perp = @(x,p)(1/175)*((-14 + 155.*C{2}(x,p) - 36*C{4}(x,p))./(1 - 0.4*C{2}(x,p)));
R.iso = @(x,p)C{2}(x,p);
    otherwise
        error('only order 1 and 2 coded');
end

%orientation
%oPara = @(tau_o,t2) (1/9).*(1+4/5.*exp(-t2./tau_o));
%oPerp = @(tau_o,t2) (1/9).*(1-2/5.*exp(-t2./tau_o));



