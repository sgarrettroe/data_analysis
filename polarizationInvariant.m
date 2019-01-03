function out = polarizationInvariant(a,b,c,d,aa,bb,cc,dd)

P = [(a'*b)*(c'*d);
    (a'*c)*(b'*d);
    (a'*d)*(b'*c)];

M = [4 -1 -1;
    -1 4 -1;
    -1 -1 4];

D = [(aa'*bb)*(cc'*dd);
    (aa'*cc)*(bb'*dd);
    (aa'*dd)*(bb'*cc)];

out = 1/30*P'*M*D;
