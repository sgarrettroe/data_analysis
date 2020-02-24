function varargout = filterExcitons(w_laser,BW,E,i_thermal,varargin)
% poor man's way to deal with finite excitation bandwidth. Just throw out
% states from the 1 and 2 exciton manifolds if they are outside the laser
% bandwidth

b = BW/2;
n = length(varargin);
for ii = 1:n
    ll = ii*w_laser - ii*b; %lower limit
    ul = ii*w_laser + ii*b; %upper limit
    ind = varargin{ii};
    i_exciton_energy = E(ind) - E(i_thermal);
    varargout{ii} = ind(i_exciton_energy>=ll & i_exciton_energy<=ul);
end
