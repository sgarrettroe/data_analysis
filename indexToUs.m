function us = indexToUs(index,lmodes)
%go from an index to a set of quantum numbers (1-based)

nstates = [lmodes.nstates];
temp = fliplr(nstates);
temp = [1 temp(1:end-1)];
factors = cumprod(temp);

n = length(lmodes);
m = size(index,1); %num rows of index

us = zeros(m,n);
    for i = 1:n
        us(:,i) = floor(index./factors(end-i+1))+1;
        index = index - ((us(:,i)-1)*factors(end-i+1));
    end
us(:,end) = us(:,end)-1;

%usToIndex(us,lmodes)
