function [O] = wobblingO(C)


O.para = @(p,t)(1/9).*(1+4/5.*C{2}(p,t));
O.perp = @(p,t)(1/9).*(1-2/5.*C{2}(p,t));