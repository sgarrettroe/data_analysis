function [out,f]=peakpos(x,y)
degree = 2;
f=polyfit(x,y,degree);
out = -f(2)/2/f(1);
