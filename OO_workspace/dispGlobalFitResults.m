function dispGlobalFitResults(gfstruct)
pfit = gfstruct.pfit;
pnames = gfstruct.pnames;

fprintf(1,'%s %s\n','fitting method: ',gfstruct.fitMethod);
fprintf(1,'%s %s\n','damping: ',gfstruct.damping);
for ii = 1:length(pfit)
    fprintf(1,'%20s\t%12f\n',pnames{ii},pfit(ii));
end