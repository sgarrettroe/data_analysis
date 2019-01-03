function timeSignalMatrix = highPassBackend(PP,t0_bin)

PP_ = PP(:,t0_bin:end); % take the elements at t >= 0
PP_var = bsxfun(@minus, PP_, mean(PP_,2)); % subtract the mean to avoid a spike in the FT at 0 frequency
R = sgrsfft(PP_var,[],2); % have to FT along the correct dimension
% excise the low frequency components
R(:,1:100) = 0;
R(:,end-100:end) = 0;

% inverse FFT and take the real part, because we added small imaginary
% components by excising the low frequency components as we did above
timeSignalMatrix = real(sgrsifft(R,[],2));