function vs = indexToVs(index,lmodes)
%go from an index to a set of quantum numbers (1-based)

us = indexToUs(index,lmodes);
vs = us-1;
