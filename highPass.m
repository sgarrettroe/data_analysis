function s2 = highPass(s)

t2 = [s.t2];
t0_bin = [s.t0_bin];

s2 = s;
for ii = 1:length(t2)
    PP = s(ii).PP;
    s2(ii).PP = highPassBackend(PP,t0_bin(ii));
end