function pmodes = evaluateAnharmonicity2(pmodes,tod,fod,varargin)

active_modes = 1:length(pmodes);
if nargin == 4
    active_modes = varargin{1};
end

F = sparse(pmodes.NSTATES,pmodes.NSTATES);

order = 3;
n_tod = size(tod,1);
for ii = 1:n_tod
    inds = tod(ii,1:order);
    val = tod(ii,order+1);
    
    %skip 
    if ~all(ismember(inds,active_modes))
        continue
    end

    % rewrite mode indicies in terms of the lmodes/scfmodes array index
    tmp = zeros(1,length(inds));
    for jj=1:length(inds)
        tmp(jj) = find(inds(jj)==active_modes);
    end
    inds = tmp;

    F_i = calcNextTerm(inds,val,pmodes);
    
    F = F + F_i;
end


order = 4;
n_fod = size(fod,1);
for ii = 1:n_fod
    inds = fod(ii,1:order);
    val = fod(ii,order+1);
    %skip 
    if ~all(ismember(inds,active_modes))
        continue
    end
    
    % rewrite mode indicies in terms of the lmodes/scfmodes array index
    tmp = zeros(1,length(inds));
    for jj=1:length(inds)
        tmp(jj) = find(inds(jj)==active_modes);
    end
    inds = tmp;

    F_i = calcNextTerm(inds,val,pmodes);
    
    F = F + F_i;
    
end

% add single mode (harmonic) hamiltonian and anh hamiltonian F
pmodes.H_ = pmodes.H + F;


end

function n = nperms(n_syms,num_repeats)
n = factorial(n_syms)/prod(factorial(num_repeats));
end

function [n_syms,n_unique,unique_syms,num_repeats] = countSymbols(in)
n_syms = length(in);

[unique_syms,~,ib]=unique(in);
n_unique = length(unique_syms);

num_repeats = zeros(1,length(unique_syms));
for ii = 1:length(ib);
    num_repeats(ib(ii)) = num_repeats(ib(ii))+1;
end

end


function [unique_syms_out,num_repeats_out] = sortSymbols(n_unique,unique_syms,num_repeats,i_mode)
% sort the list of modes in terms of the active mode first and then the
% rest, in parallel sort the modes that are getting the mean-field
% treatment
num_repeats_out = zeros(1,n_unique);
unique_syms_out = zeros(1,n_unique);
ind_mode = (unique_syms==i_mode);

num_repeats_out(1) = num_repeats(ind_mode);
unique_syms_out(1) = unique_syms(ind_mode);

num_repeats_out(2:end) = num_repeats(~ind_mode);
unique_syms_out(2:end) = unique_syms(~ind_mode);

end

function  F_i = calcNextTerm(inds,val,pmodes)
[n_syms,n_unique,unique_syms,num_repeats] = countSymbols(inds);
%n_perms = nperms(n_syms,n_repeats);
%prefactor = 1/factorial(n_syms);
prefactor = 1/prod(factorial(num_repeats));
[unique_syms,n] = sortSymbols(n_unique,unique_syms,num_repeats,inds(1));

OP = speye(pmodes.NSTATES,pmodes.NSTATES);

for i_mode = 1:n_unique
    
fname = sprintf('Q_POW_%i%i',n(i_mode),unique_syms(i_mode));
OP = OP*pmodes.(fname);

end
%F_i = prefactor*n_perms*val*OP;
F_i = prefactor*val*OP;
end