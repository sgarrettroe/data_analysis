function maxMatrix = CLSMaxMatrix(peakFit)

[m,n] = size(peakFit);
maxMatrix = zeros(m,n,2);
for ii = 1:m
    for jj = 1:n
         maxMatrix(ii,jj,1) = peakFit(ii,jj).fitresult.center;
         dummy = confint(peakFit(ii,jj).fitresult);
         dummy_coeff_names = coeffnames(peakFit(ii,jj).fitresult);
         centerFreqIndex = ~cellfun('isempty',strfind(dummy_coeff_names, 'center'));
         err_b = dummy(:,centerFreqIndex);
         maxMatrix(ii,jj,2) = (err_b(2)-err_b(1))/2;
    end
end