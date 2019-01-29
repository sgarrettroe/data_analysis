function viewCLS(CLSstruct,ii)
if isscalar(ii)
%     figure,clf
    hold on
    errorbar(CLSstruct(ii).w1,CLSstruct(ii).center,CLSstruct(ii).center_std,'r.')
    plot(CLSstruct(ii).w1,CLSstruct(ii).fitresult(CLSstruct(ii).w1))
    hold off
else
    warning('Non-scalar input for scan number')
    return
end