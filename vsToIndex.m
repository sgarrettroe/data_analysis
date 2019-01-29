function index = vsToIndex(vs,lmodes)
%convert input quantum numbers of states (0-based) and convert that to the number
%of the product mode basis vector.

us = vs+1;
index = usToIndex(us,lmodes);

