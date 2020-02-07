function s = stats(varargin)
%STATS.M calculate the moments and cumulants 1-3 of a distribution
% call as
%
% s = stats(x,p) where x, and p are arrays of an axis and a
% distribution. The output is a structure of moments and cumulants
%
% s = stats(@p,llim,ulim)
%

%fun = @trapz;
fun = @sum;

if isa(varargin{1},'function_handle')
  p = varargin{1};
  llim = varargin{2};
  ulim = varargin{3};
  
  parse_options(varargin{4:end});
  
  norm = quad(p,llim,ulim);
  m1 = quad(@(x) x.*p(x),llim,ulim);
  m2 = quad(@(x) x.^2.*p(x),llim,ulim);
  m3 = quad(@(x) x.^3.*p(x),llim,ulim);
else
  x = varargin{1};
  p = varargin{2};
  
  parse_options(varargin{3:end});

  dx = x(2)-x(1);
  norm = trapz(p)*dx;
  
  %moments of x
  m1 = feval(fun,x.*p)*dx/norm;
  m2 = feval(fun,x.^2.*p)*dx/norm;
  m3 = feval(fun,x.^3.*p)*dx/norm;
%   m2 = trapz(x.^2.*p)*dx/norm;
%   m3 = trapz(x.^3.*p)*dx/norm;
  
end

%cumulants or "central moments"
% these are variance, skewness
c1 = m1;
c2 = m2 - m1^2;
c3 = m3 - 3*m1*c2 - m1^3;

s.norm = norm;
s.m1 = m1;
s.m2 = m2;
s.m3 = m3;
s.c1 = c1;
s.c2 = c2;
s.c3 = c3;

  function parse_options(varargin)
    while length(varargin)>=2
      arg = varargin{1};
      val = varargin{2};
      switch lower(arg)
        case {'fun'}
          %fun = eval(['@' val]);
          fun = val;
        otherwise
          error(['unknown option ',arg])
      end
      varargin = varargin(3:end);
    end
  end %end nested function
end %end main function
