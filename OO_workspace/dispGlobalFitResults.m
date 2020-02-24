function dispGlobalFitResults(gfstruct,varargin)
% dispGlobalFitResults(gfstruct) outputs the results of global fitting to
% the command window, given an input structure 'GFSTRUCT' that contains the
% fitting method, damping, pnames, and pfit.
%
% If using the standard method for your global fit structure, this would be
% found in fitStruct(ii,jj,...).RESULTS(kk).
fid = 1;

while length(varargin)>=2 %using a named pair
    arg = varargin{1};
    val = varargin{2};
    switch lower(arg)
        case 'fid'
            fid = val;
        otherwise
            warning(['dispGlobalFitResults: unknown option ',arg])
    end
    varargin = varargin(3:end);
end

pfit = gfstruct.pfit;
pnames = gfstruct.pnames;

fprintf(fid,'%s %s\n','fitting method: ',gfstruct.fitMethod);
fprintf(fid,'%s %s\n','damping: ',gfstruct.damping);
for ii = 1:length(pfit)
    fprintf(fid,'%20s\t%12f\n',pnames{ii},pfit(ii));
end