function index = usToIndex(us,lmodes)
%convert input numbers of states (1-based) and convert that to the number
%of the product mode basis vector.

nstates = [lmodes.nstates];
temp = fliplr(nstates);
temp = [1 temp(1:end-1)];
factors = cumprod(temp);

us = flipud(us(:));
index = factors*(us-1)+1;