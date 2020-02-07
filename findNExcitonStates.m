function varargout = findNExcitonStates(PSIi,C0,sig_figs)
% find the indicies into the wavefunction for the 1:N exciton states. The
% results are returned in a cell array (because they are all different
% sizes). The creation operator is used to excite the initial wavefunction,
% PSIi. Probably best to do all this in the eigenstate basis. For example
% after:
%
% [V,E]=eig(H_,'vector'); [E,ordering] = sort(E); V = V(:,ordering);
% %eigenvectors in input basis VV = eye(size(V)); %eigenvectors in
% eigenstate basis PSIi = VV(:,1); %take first eigenstate for the time
% being C0 = V'*C*V;
%
% then
%
% [ind1 ind2 ... indn ] = findNExcitonStates(PSIi,C0,2);
% 
% will find up to the nth exciton states using a 2 digits after the decimal
% (hundredths) resolution. (so only states with more than 1% amplitude will
% be included.) Preliminarily it looks like 1 might be better...
%

varargout = cell(1,nargout);
for ii = 1:nargout
    varargout{ii} = find(abs(round(C0^ii*PSIi,sig_figs))>0);
end
