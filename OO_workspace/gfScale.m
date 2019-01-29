function [scale] = gfScale(dataMatrix)
[m,n,~] = size(dataMatrix);
scale = repmat(abs(min(min(dataMatrix,[],2),[],1)),m,n);