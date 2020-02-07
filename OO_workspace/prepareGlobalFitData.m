function [dataMatrix,sdMatrix,weightsMatrix] = prepareGlobalFitData(croppedData)
% [dataMatrix,sdMatrix,weightsMatrix] = prepareGlobalFitData(croppedData)
%
% prepareGlobalFitData takes an input of a cropped data structure and
% returns a M x N x P matrix that can be used for global fitting of data.

[m,n] = size(croppedData(1).R);
p = length(croppedData);

% Replaced a couple of FOR loops that used to define our data matrices by
% looping over CROPPEDDATA, and find the MINIMUM of each matrix R in
% CROPPEDDATA, and then NORMALIZE the data by the scale factor.
%
% Increased the speed of the code by a factor of ~37 by doing so.
%
% Sorry for the MatLab gibberish though.


dataMatrix = reshape([croppedData(:).R],m,n,p);

% 2016-06-20 We now handle defining the scale directly in any global
% fitting method that needs it
% scale = repmat(abs(min(min(dataMatrix,[],2),[],1)),m,n);
% dataMatrix = dataMatrix./scale;
if isfield(croppedData,'sd')
    sdMatrix = reshape([croppedData(:).sd],m,n,p);
    varMatrix = sdMatrix.^2;
    weightsMatrix = 1./varMatrix;
end