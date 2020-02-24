function pmodes = evaluateAnharmonicity(pmodes,tod,fod,varargin)
% evaluateAnharmonicity.m adds third and fourth-order anharmonic constants
% to the vibrational hamiltonian in the product mode basis.
%
% pmodes = evaluateAnharmonicity(pmodes,tod,fod,active_modes) takes the
%      modes in a product basis, pmodes, the third-order derivatives, tod,
%      and fourth-order derivatives, fod, and the list of active_modes
%
% pmodes = evaluateAnharmonicity(pmodes,tod,fod) assumes that all modes are
% active.

active_modes = 1:length(pmodes);
if nargin == 4
    active_modes = varargin{1};
end

n_tod = size(tod,1);
for ii = 1:n_tod
    inds = tod(ii,1:3);
    val = tod(ii,4);
    %skip 
    if ~all(ismember(inds,active_modes))
        continue
    end
    ind1 = find(inds(1)==active_modes);
    ind2 = find(inds(2)==active_modes);
    ind3 = find(inds(3)==active_modes);
    P = unique(perms(inds),'rows');
    n_perms = size(P,1);
    fname1 = sprintf('Q%i',ind1);
    fname2 = sprintf('Q%i',ind2);
    fname3 = sprintf('Q%i',ind3);
    OP1 = pmodes.(fname1);
    OP2 = pmodes.(fname2);
    OP3 = pmodes.(fname3);
 
    pmodes.H_ = pmodes.H_ + 1/6*n_perms*val*OP1*OP2*OP3;
end

n_fod = size(fod,1);
for ii = 1:n_fod
    inds = fod(ii,1:4);
    val = fod(ii,5);
    %skip 
    if ~all(ismember(inds,active_modes))
        continue
    end
    ind1 = find(inds(1)==active_modes);
    ind2 = find(inds(2)==active_modes);
    ind3 = find(inds(3)==active_modes);
    ind4 = find(inds(4)==active_modes);
    P = unique(perms(inds),'rows');
    n_perms = size(P,1);
    fname1 = sprintf('Q%i',ind1);
    fname2 = sprintf('Q%i',ind2);
    fname3 = sprintf('Q%i',ind3);
    fname4 = sprintf('Q%i',ind4);
    OP1 = pmodes.(fname1);
    OP2 = pmodes.(fname2);
    OP3 = pmodes.(fname3);
    OP4 = pmodes.(fname4);
 
    pmodes.H_ = pmodes.H_ + 1/24*n_perms*val*OP1*OP2*OP3*OP4;
end

%add the harmonic hamiltonian back in also
pmodes.H_ = pmodes.H + pmodes.H_;
