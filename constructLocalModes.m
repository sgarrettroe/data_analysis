function lmodes = constructLocalModes(options)
% constructLocalModes takes an option structure array in and outputs a structure
% array with operators and basis states for each local mode

hbar = 1; %might want to change later?

% set up structure for output
lmodes = struct('nstates',[],...
    'identity',[],...
    'u',[],...
    'v',[],...
    'a',[],...
    'c',[],...
    'q',[],...
    'p',[],...
    'mux',[],...
    'muy',[],...
    'muz',[],...
    'h',[],...
    'h_',[]);

for ii = 1:length(options)
    nstates = options(ii).nstates;
    m = options(ii).m;
    w = options(ii).w;%*2*pi*2.9979e10/1e-12; %convert cm-1 to rad/ps
    dmu_dq = options(ii).dmu_dq(:); %should be a column vector of x, y, z components of the dipole
    
    %creation and annihil ops have amplitudes that scale like sqrt(n)
    prefactors = sqrt(1:nstates-1);
    
    %identity ops for each subspace
    identity = eye(nstates);
    
    %
    % mode ii
    %
    
    % basis for mode 1 call them u
    u = @(ii) identity(:,ii);
    
    % projectors for mode 1
    proj = @(ii) u(ii)'*u(ii);
    
    % shortcut for v=0 type notation
    v = @(ii) u(ii+1);
    
    % annihilation operator
    a = diag(prefactors,1);
    % creation operator
    c = diag(prefactors,-1);
    
    % position operator
    q = sqrt(hbar/(2*m*w)).*(a + c);
    
    % momentum operator
    p = -1i*sqrt(hbar*m*w/2).*(a - c);
    
    % hamiltonian
    h = hbar*w*(c*a + 1/2*identity);
    h_ = p*p/(2*m) + m*w^2*q*q/2;
    
    % dipole operator
    mux = dmu_dq(1).*(a + c);
    muy = dmu_dq(2).*(a + c);
    muz = dmu_dq(3).*(a + c);
    
    lmodes(ii).nstates=nstates;
    lmodes(ii).identity=identity;
    lmodes(ii).u=u;
    lmodes(ii).v=v;
    lmodes(ii).a=a;
    lmodes(ii).c=c;
    lmodes(ii).q=q;
    lmodes(ii).p=p;
    lmodes(ii).mux=mux;
    lmodes(ii).muy=muy;
    lmodes(ii).muz=muz;
    lmodes(ii).h=h;
    lmodes(ii).h_ = h_;
    
end