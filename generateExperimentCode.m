function out = generateExperimentCode(varargin)
% generateExperimentCode Create a random code for experiment tracking.
%   C = generateExperimentCode provides a random alphanumeric code using
%   the default version.
%
%   C = generateExperimentCode('version',v) generates the code in version
%   v.
%
%   Version 1: four digit uppercase alpha plus four digit numeric code,
%   e.g. ABCD1234.

% parse input args
default_version = 1;
MAX_VERSION = 1;

p = inputParser;
validVersionNumber = @(x) isnumeric(x) & x<=MAX_VERSION;
addParameter(p,'version',default_version,validVersionNumber);
parse(p,varargin{:});

version  = p.Results.version;

switch version
    case 1
        f = @codeVersion1;
    otherwise
        error('unknown version number %f',version)
end
out = f;
end

function out = codeVersion1
% codeVersion1 generate codes like ABCD1234
n_letters = 4;
n_numbers = 4;
letter_list = 'A':'Z';
number_list = '0':'9';
out = [letter_list(randi(numel(letter_list),[1,n_letters])),...
    number_list(randi(numel(number_list),[1,n_numbers]))];
end