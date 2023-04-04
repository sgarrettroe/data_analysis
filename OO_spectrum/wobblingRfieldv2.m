function [R] = wobblingRfieldv2(C,stark_effect_order)
% wobbling in a cone model of Kramer Nishda Fayer JCP 2016
% Field accomodations are from Hoffman2021

S{1} = @(p) coth(p.BetaV)-1/p.BetaV;
S{2} = @(p) (3*p.BetaV*coth(p.BetaV-p.BetaV^2-3))/p.BetaV^2;
S{3} = @(p) ((p.BetaV^2+15)*(p.BetaV)*coth(p.BetaV-6*p.BetaV^2-15))/p.BetaV^3; %Actually unnecessary unless you want to show the t->inf points
 
a = @(p) 1/3 + 2/3*(S{2}(p));
b = @(p) (S{1}(p))^2;

switch stark_effect_order
    case 1 % from Kramer2015a
R.paraInt = @(x,p)(3/25).*((11.*C{1}(x,p) + 4.*C{3}(x,p))./(1 + 0.8.*C{2}(x,p)));
R.perpInt = @(x,p)(3/25).*(( 7.*C{1}(x,p) - 2.*C{3}(x,p))./(1 - 0.4.*C{2}(x,p)));
R.isoInt  = @(x,p)C{1}(x,p);
    case 2 % from 
R.paraInt = @(x,p)(1/175).*(( 28 + 215.*C{2}(x,p) + 72.*C{4}(x,p))./(1 + 0.8.*C{2}(x,p)));
R.perpInt = @(x,p)(1/175).*((-14 + 155.*C{2}(x,p) - 36.*C{4}(x,p))./(1 - 0.4.*C{2}(x,p)));
R.isoInt  = @(x,p)C{2}(x,p);
    otherwise
        error('only order 1 and 2 coded');       
end

R.para = @(x,p) ((1-b(p))*R.paraInt(x,p)+(b(p)-b(p)/a(p)).*C{1}(x,p)-b(p)+(b(p)^2)/a(p))/((1-b(p))*(1-b(p)/a(p)));
R.perp = @(x,p) ((1-b(p))*R.perpInt(x,p)+(b(p)-b(p)/a(p)).*C{1}(x,p)-b(p)+(b(p)^2)/a(p))/((1-b(p))*(1-b(p)/a(p)));
R.iso  = @(x,p) ((1-b(p))*R.isoInt(x,p) +(b(p)-b(p)/a(p)).*C{1}(x,p)-b(p)+(b(p)^2)/a(p))/((1-b(p))*(1-b(p)/a(p)));
end 
