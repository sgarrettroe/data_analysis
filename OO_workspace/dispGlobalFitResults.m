function dispGlobalFitResults(gfstruct)
% dispGlobalFitResults(gfstruct) outputs the results of global fitting to
% the command window, given an input structure 'GFSTRUCT' that contains the
% fitting method, damping, pnames, and pfit.
%
% If using the standard method for your global fit structure, this would be
% found in fitStruct(ii,jj,...).RESULTS(kk).
pfit = gfstruct.pfit;
pnames = gfstruct.pnames;

fprintf(1,'%s %s\n','fitting method: ',gfstruct.fitMethod);
fprintf(1,'%s %s\n','damping: ',gfstruct.damping);
for ii = 1:length(pfit)
    fprintf(1,'%20s\t%12f\n',pnames{ii},pfit(ii));
end