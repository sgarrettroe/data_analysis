function [dataobj] = cropData(data,range1,range3)
% CROPDATA returns a cropped section of 2D-IR data.
%
%   CROPPEDDATAOBJ = cropData(DATASTRUCT,w1_range,w3_range)
%   creates a 1xn data structure array, DO, the
%   contains a cropped version of the w1, w3, and response values from
%   2D-IR data (based on the input ranges specified), as well as the t2
%   times.
%
%   --DATASTRUCT must be a a 1xn 2D-IR structure array with the fields
%   data(n).w1,data(n).w3, data(n).R (the response), data(n).t2 (population
%   times)
%
%   --w1_range should be the w1 range in cm-1 [w1(1) w1(end)]
%   --w3_range should be the w3 range in cm-1 [w3(1) w3(end)]
%         
%   This function was designed to allow (1) the ability to omit indices
%   from the fitting functions we use by cropping data to the desired size
%   initially (or in the course of the function), and (2) to reduce the
%   size of the data object being passed through successive methods in
%   fitting, since the resulting data object is approximately 0.1% of the
%   size of the original data array for typical data.


dataobj = struct('w1',[],'w3',[],'R',[]);
for ii = 1:length(data)
    if isempty(data(ii).R)
        continue
    end
    ind1 = find(data(ii).w1>range1(1) & data(ii).w1<range1(2));
    ind3 = find(data(ii).w3>range3(1) & data(ii).w3<range3(2));
    dataobj(ii).w1 = data(ii).w1(ind1);
    dataobj(ii).w3 = data(ii).w3(ind3);
    dataobj(ii).R = data(ii).R(ind3,ind1);
    dataobj(ii).t2 = data(ii).t2;
end