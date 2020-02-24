function out = interp2D(in,zp_factor)
% 'in' needs to be a data structure that we've generated from our 2D.
% 'zp_factor' is the scalar you're going to multiply the n_zp (number of
% zero padding) by. Increasing this scalar will decrease the w1 spacing in
% the output compared with the input.
out = in;
for ii = 1:length(in)
    in(ii).zeropad = zp_factor*in(ii).zeropad;
    % a placeholder to later restore our w3 calibration
    w3 = in(ii).w3; 
    
    % the heart of the calcuation
    s = absorptive2dPP(in(ii),'zeropad',in(ii).zeropad);
    out(ii) = s;
    
    %fix the sign on our output
    out(ii).R = -out(ii).R;
    
    % redefine w3 to match our previous calibration, which we lost when we
    % reran the ABSORPTIVE2DPP function
    out(ii).w3 = w3; 
end