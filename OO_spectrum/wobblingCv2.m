function [C,S,tau_eff] = wobblingCv2
% wobbling in a cone model of Kramer Nishda Fayer JCP 2016

%stark_effect_order = varargin{1};

%tau_eff = in.tau;
%theta_deg = in.theta_deg;

x0 = @(theta_deg) cos(theta_deg/180*pi);

S{1} = @(theta_deg) (1+x0(theta_deg))/2;
S{2} = @(theta_deg)  x0(theta_deg)*(1+x0(theta_deg))/2;
S{3} = @(theta_deg) (1+x0(theta_deg))*(5*x0(theta_deg)^2-1)/8;
S{4} = @(theta_deg)  x0(theta_deg)*(1+x0(theta_deg))*(7*x0(theta_deg)^2-3)/8;

D = @(tr)1/(6*tr);
theta = @(theta_deg) theta_deg*pi/180;

% small angle approx
%tau_eff{1} = @(tr,theta_deg) theta(theta_deg)^2/D(tr)*(7/24-55*theta(theta_deg)^2/1152);
%tau_eff{2} = @(tr,theta_deg) theta(theta_deg)^2/D(tr)*(7/24-35*theta(theta_deg)^2/384);
%tau_eff{3} = @(tr,theta_deg) theta(theta_deg)^2/D(tr)*(7/24-5*theta(theta_deg)^2/32);
%tau_eff{4} = @(tr,theta_deg) theta(theta_deg)^2/D(tr)*(7/24-35*theta(theta_deg)^2/144);
%sm = tau_eff;
 
tau_eff{1} = @(tr,theta_deg) (-(1+x0(theta_deg))^2/(2*(1-x0(theta_deg)))*log((1+x0(theta_deg))/2)...
    -x0(theta_deg)*(4+x0(theta_deg)-x0(theta_deg)^2)/4)...
    /(D(tr)*(1-S{1}(theta_deg)^2));

tau_eff{2} = @(tr,theta_deg) (-x0(theta_deg)^2*(1+x0(theta_deg))^2*(log((1+x0(theta_deg))/2)+(1-x0(theta_deg))/2)...
    /(2*(1-x0(theta_deg)))...
    +(1-x0(theta_deg))*(6+8*x0(theta_deg)-x0(theta_deg)^2-12*x0(theta_deg)^3-7*x0(theta_deg)^4)/24 ...
    )...
    /(D(tr)*(1-S{2}(theta_deg)^2));

tau_eff{3} = @(tr,theta_deg) (-(1+x0(theta_deg))^2*(5*x0(theta_deg)^2-1)^2*log((1+x0(theta_deg))/2)...
    /(32*(1-x0(theta_deg)))...
    +(30-9*x0(theta_deg)+154*x0(theta_deg)^2+231*x0(theta_deg)^3-370*x0(theta_deg)^4-615*x0(theta_deg)^5+10*x0(theta_deg)^6+185*x0(theta_deg)^7)...
    /384 ...
    )...
    /(D(tr)*(1-S{3}(theta_deg)^2));

%tau_eff{4} is left as an exercise for the reader... sorry.
%reader: rude
tau_eff{4} = @(tr,theta_deg) (-x0(theta_deg)^2*(1+x0(theta_deg))^2*(7*x0(theta_deg)^2-3)^2*log((1+x0(theta_deg))/2)...
    /(32*(1-x0(theta_deg)))...
    +(260+68*x0(theta_deg)-1301*x0(theta_deg)^2-1895*x0(theta_deg)^3+6231*x0(theta_deg)^4+9597*x0(theta_deg)^5-7511*x0(theta_deg)^6-13517*x0(theta_deg)^7 + ...
    497*x0(theta_deg)^8 + 3731*x0(theta_deg)^9)...
    /3840 ...
    )...
    /(D(tr)*(1-S{4}(theta_deg)^2));

C = cell(1,4);

for l = 1:4
    
    C{l} = @(t,tr,theta_deg) S{l}(theta_deg)^2 + (1-S{l}(theta_deg)^2)*exp(-t./tau_eff{l}(tr,theta_deg));
    
end


