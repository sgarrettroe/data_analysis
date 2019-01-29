function r = orientationalResponse(tau_c,order,varargin)
% orientationalResponse.m Calculate the value of the orientational
% correlation function for different orders of spectroscopy *for the all
% parallel case only*!
%
% call as r = orientationalResponse(tau_c,order,[times])
% so for 1st order
% r = orientationalResponse(tau_c,order,t)
%
% For 3rd order
% r = orientationalResponse(tau_c,order,T1,t2,T3)
%
% For 3rd order
% r = orientationalResponse(tau_c,order,T1,t2,T3,t4,T5)
%
% tau = 1/(6D), i.e. the value from pump-probe anistropy decay measurements
%   just be careful to use the same units as the time axis/axes
%
% order = [1,3,5]
%
%
% mostly based on a handout from Jan Helbing, which is largely based on the
% papers by Tokmakoff JCP 105, 1-12 (1996) and JCP 105, 13 (1996)
% paper by Helbing JCP 122, 124505 (2005)
% Zanni also has a paper Ding et al JCP 123, 094592 (2005)
D = 1/(6*tau_c);
c1 = @(t) exp(-2*D.*t);
c2 = @(t) exp(-6*D.*t);
c3 = @(t) exp(-12*D.*t);
if length(varargin)~=order
  error('Orientational response: wrong number of input time axes')
end
switch order
  case 1
    t = varargin{1};
    r = exp(-2*D.*t);
  case 3
    T1 = varargin{1};
    t2 = varargin{2};
    T3 = varargin{3};
    r = 1/9.*c1(T1).*(1 + 4/5.*c2(t2)).*c1(T3);
  case 5
    T1 = varargin{1};
    t2 = varargin{2};
    T3 = varargin{3};
    t4 = varargin{4};
    T5 = varargin{5};
    r = c1(T1).*c1(T5).*((1/27)*(1+4/5.*c2(t2)).*(1+4/5.*c2(t4)).*c1(T3)...
      +4/175*c2(t2)*c2(t4).*c3(T3));
  otherwise
    error(['orientationalResponse.m does not support order ',num2str(order),' only 1,3,5'])
end
