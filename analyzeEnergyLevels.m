function [V,E,ind] = analyzeEnergyLevels(lmodes,pmodes,options)
arguments
    lmodes (1,:) struct
    pmodes (1,1) struct
    options.roptions (1,1) struct = struct()
    options.BW (1,1) double = realmax
    options.ind (1,:) double = []  % candidate for deletion
    options.method char {mustBeMember(options.method,{'vci','vpt2'})} = 'vci'
    options.n_exciton_sig_figs (1,1) double = 1
    options.n_sparse_states (1,1) double = min(pmodes.NSTATES-2,200)
    options.w_laser (1,1) double = 2000
end

% named/keyval pairs first
BW = options.BW;
ind = options.ind;  % candidate for deletion
n_exciton_sig_figs = options.n_exciton_sig_figs;
n_sparse_states = options.n_sparse_states;
w_laser = options.w_laser;

% overwrite with roptions if provided
if isfield(options.roptions, 'ind')  % candidate for deletion
    ind = options.roptions.ind;  
end
if isfield(options.roptions, 'BW')
    BW = options.roptions.BW;
end
if isfield(options.roptions, 'n_exciton_sig_figs')
    n_exciton_sig_figs = options.roptions.n_exciton_sig_figs;
end
if isfield(options.roptions, 'n_sparse_states')
    n_sparse_states = options.roptions.n_sparse_states;
end
if isfield(options.roptions, 'w_laser')
    w_laser = options.roptions.w_laser;
end

if ~isempty(ind)
    warning('SGRLab:Anharm:analyzeEnergyLevels',...
        'option "ind" is deprecated. States are chosen automatically.')
end

% unpack needed matrices
C = pmodes.C;
H = pmodes.H;
H_ = pmodes.H_;
MUX = pmodes.MUX;
MUY = pmodes.MUY;
MUZ = pmodes.MUZ;

switch lower(options.method)
    case 'vci'
        [V, E, original_energies, original_vecs] = ...
            anharm_method_vci(H, H_, n_sparse_states);
    case 'vpt2'
        if ~isfield(pmodes, 'H__')
            error('SGRLab:vpt2ProductModes:missingValue', ...
                ['Input pmodes must have field H__. ', ...
                'Try running\n', ...
                '>> pmodes = evaluateAnharmonicity2(pmodes,tod,fod);']);
        end
        H__ = pmodes.H__;
        [V, E, original_energies, original_vecs] = ...
            anharm_method_vpt2(H, H__, n_sparse_states);
end

% VV is eigenvectors in eigenstate basis
if issparse(H_)
    VV = speye(length(E),length(E)); 
else
    VV = eye(size(V));
end

if issparse(original_energies)
    original_energies = full(original_energies);
end

new_energies = full(E);

ind_0ex = 1; %take only the first state as the ground state
C0 = V'*C*V;
PSIi = VV(:,ind_0ex); %take first eigenstate for the time being

%find all one and two exciton states
[ind_1ex, ind_2ex] =  findNExcitonStates(PSIi,C0,n_exciton_sig_figs);

%keep only the ones in the laser bandwidth
[ind_1ex, ind_2ex] = filterExcitons(w_laser,BW,E,1,ind_1ex,ind_2ex);

ind = {ind_0ex;ind_1ex;ind_2ex};

% calculate energy differences relative ot the ground state
original_gaps = original_energies - original_energies(ind_0ex);
new_gaps = new_energies - new_energies(ind_0ex);

%for ii = 1:min(6,length(V))
%print energy of all states in gs, 1 ex, and 2
disp('States by energy');
iis = unique(vertcat(ind{:})'); %the ' is important to make iis a row vec
for ii = iis
    fprintf(1,'i = %3d\tE0 = %-8.1f E = %-8.1f gap0 = %-8.1f gap = %-8.1f \n',...
        ii,original_energies(ii),new_energies(ii),original_gaps(ii),new_gaps(ii));
end
fprintf(1,'\n');

%
iis = unique(vertcat(ind{:})');
for ii = iis
    fprintf(1,'i = %3d\tenergy above g.s. = %8.1f\n', ii, new_gaps(ii));
    displayCoeffMatrix(V(:,ii),lmodes);
end


for ii = iis
    eigenstate_of_interest = ii;
    vec = V(:,eigenstate_of_interest);
    vec = conj(vec).*vec;
    [~,ind] = sort(vec,'descend');
    vec = vec(ind);
    ind = ind(vec>=0.1);
    vec = vec(vec>=0.1);
    
    vs = indexToVs(ind,lmodes);
    
    %set up some output
    nmodes = length(lmodes);
    
    formatstring = '(';
    for jj = 1:nmodes
        formatstring = strcat(formatstring,'%3d,');
    end
    formatstring(end) = ')';
    
    fprintf(1,'EIGENSTATE %3d\n',eigenstate_of_interest);
    fprintf(1,'index\tcontribution\tquantum nums\n');
    for jj = 1:length(ind)
        fprintf(1,strcat('%d\t%f\t',formatstring,'\n'),ind(jj),vec(jj));
        for kk = 1:nmodes
            fprintf(1,'%3d,',vs(jj,kk));
        end
        %\b is a backspace... 
        fprintf(1,'\b)\n');
    end
    fprintf(1,'\n');
end
fprintf(1,'\n');

% transitions starting from ground state
disp('ZERO TO ONE EXCITON TRANSITIONS');
fprintf('i\tj\to_gap\tn_gap\tu_orig\tu_mixed\n');
iis = ind_0ex';
jjs = ind_1ex';
omu = zeros(1,length(jjs));
nmu = zeros(1,length(jjs));
for ii = iis
    for jj = jjs
        % original states
        PSIi = original_vecs(:,ii);
        PSIf = original_vecs(:,jj);
        
        temp = [PSIf'*MUX*PSIi;...
            PSIf'*MUY*PSIi;...
            PSIf'*MUZ*PSIi];
        
        omu(jj) = sqrt(temp'*temp);
        
        %mixed states
        PSIi = V(:,ii);
        PSIf = V(:,jj);
        
        temp = [PSIf'*MUX*PSIi;...
            PSIf'*MUY*PSIi;...
            PSIf'*MUZ*PSIi];

        nmu(jj) = sqrt(temp'*temp);

        fprintf(1,'%-8d%-8d%-8.1f%-8.1f%-8.1f%-8.1f\n',ii,jj,...
            original_gaps(jj),new_gaps(jj),omu(jj),nmu(jj));
        
    end
end

disp('')
disp('ONE TO TWO EXCITON TRANSITIONS');
fprintf('i\tj\to_gap\tn_gap\tu_orig\tu_mixed\n');
iis = ind_1ex';
jjs = ind_2ex';
omu = zeros(1,length(jjs));
nmu = zeros(1,length(jjs));
for ii = iis
    for jj = jjs
        % original states
        PSIi = original_vecs(:,ii);
        PSIf = original_vecs(:,jj);
        
        temp = [PSIf'*MUX*PSIi;...
            PSIf'*MUY*PSIi;...
            PSIf'*MUZ*PSIi];
        
        omu(jj) = sqrt(temp'*temp);
        
        %mixed states
        PSIi = V(:,ii);
        PSIf = V(:,jj);
        
        temp = [PSIf'*MUX*PSIi;...
            PSIf'*MUY*PSIi;...
            PSIf'*MUZ*PSIi];

        nmu(jj) = sqrt(temp'*temp);

        fprintf(1,'%-8d%-8d%-8.1f%-8.1f%-8.1f%-8.1f\n',ii,jj,...
            original_energies(jj)-original_energies(ii),...
            new_energies(jj)-new_energies(ii),omu(jj),nmu(jj));
        
    end
end

end

function [V, E, original_energies, original_vecs] = anharm_method_vci(H, H_, n_sparse_states)
%anharm_method_vci Diagonalize hamiltonian in product mode basis.

if issparse(H_)
    [original_vecs,original_energies] = eigs(H,n_sparse_states,'SM');
    original_energies = diag(original_energies);
    [original_energies,ordering] = sort(original_energies);
    original_vecs = original_vecs(:,ordering);

    % calculate the eigenvectors of the coupled system
    [V,E]=eigs(H_,n_sparse_states,'SM');
    E = diag(E);
    [E,ordering] = sort(E);
    V = V(:,ordering); %eigenvectors in input basis
else
    [original_vecs,original_energies] = eig(H,'vector');
    [original_energies,ordering] = sort(original_energies);
    original_vecs = original_vecs(:,ordering);
    % calculate the eigenvectors of the coupled system
    [V,E]=eig(H_,'vector');
    [E,ordering] = sort(E);
    V = V(:,ordering); %eigenvectors in input basis
end

end

function [V, E, original_energies, original_vecs] = anharm_method_vpt2(H, H__, n_sparse_states)
%anharm_method_vpt2 Second-order vibrational perturbation theory

[NSTATES, ~] = size(H);

% zeroth-order energy and vectors
E_0 = diag(H);
V_0 = speye(NSTATES);

% first-order energy and vectors
E_1 = zeros(NSTATES,1);
V_1 = zeros(NSTATES, NSTATES);

% second-order energy and vectors
E_2 = zeros(NSTATES,1);
V_2 = zeros(NSTATES, NSTATES);

for nn = 1:NSTATES

    H__nn = H__(nn,nn);
    E_1(nn) = H__nn;

    for kk = 1:NSTATES

        if nn == kk, continue, end

        H__kn = H__(kk,nn);
        Delta_E_nk = E_0(nn) - E_0(kk);

        % first-order correction to the wavefunction
        V_1(kk,nn) = V_1(kk,nn) + H__kn / Delta_E_nk;

        % second-order correction to the energy
        E_2(nn) = E_2(nn) + H__kn * conj(H__kn) / Delta_E_nk;

        % second-order correction to the wavefunction
        V_2(kk,nn) = V_2(kk,nn) - H__kn * H__nn / Delta_E_nk^2;
        V_2(nn,nn) = V_2(nn,nn) - 0.5 * H__kn * conj(H__kn) / Delta_E_nk^2;

        for ll = 1:NSTATES
            if nn == ll, continue, end

            Delta_E_nl = E_0(nn) - E_0(ll);
            H__kl = H__(kk,ll);
            H__ln = H__(ll,nn);

            V_2(kk,nn) = V_2(kk,nn) ...
                + H__kl * H__ln / (Delta_E_nk * Delta_E_nl);
        end

    end

end

% calc total energy and corrected wavefunction
V = V_0 + V_1 + V_2;
E = E_0 + E_1 + E_2;

% put everything in energy order
[E,ordering] = sort(E);
V = V(:,ordering); %eigenvectors in input basis
[original_energies, ordering] = sort(E_0);
original_vecs = V_0(:, ordering);
end
