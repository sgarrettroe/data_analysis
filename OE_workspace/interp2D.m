function out = interp2D(in,zp_factor)
out = in;
for ii = 1:length(in)
    in(ii).zeropad = zp_factor*in(ii).zeropad;
    s = absorptive2dPP(in(ii),'zeropad',in(ii).zeropad);
    out(ii) = s;
end