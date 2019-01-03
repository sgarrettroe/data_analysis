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
    'q_pow_1',[],...
    'q_pow_2',[],...
    'q_pow_3',[],...
    'q_pow_4',[],...
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
    %q = sqrt(hbar/(2*m*w)).*(a + c);
    q = sqrt(1/2).*(a + c);
    
    % momentum operator
    %p = -1i*sqrt(hbar*m*w/2).*(a - c);
    p = -1i*sqrt(1/2).*(a - c);
    
    % hamiltonian
    h = hbar*w*(c*a + 1/2*identity);
    %h_ = p*p/(2*m) + m*w^2*q*q/2;
    h_ = hbar*w*(p*p/2 + q*q/2);
    
    % dipole operator
    mux = dmu_dq(1).*(a + c);
    muy = dmu_dq(2).*(a + c);
    muz = dmu_dq(3).*(a + c);
    
    out = matrixElements_q1_q2_q3_q4(nstates);
    
    lmodes(ii).nstates=nstates;
    lmodes(ii).identity=identity;
    lmodes(ii).u=u;
    lmodes(ii).v=v;
    lmodes(ii).a=a;
    lmodes(ii).c=c;
    lmodes(ii).q=q;
    lmodes(ii).p=p;
    lmodes(ii).q_pow_1 = out.x1;
    lmodes(ii).q_pow_2 = out.x2;
    lmodes(ii).q_pow_3 = out.x3;
    lmodes(ii).q_pow_4 = out.x4;
    lmodes(ii).mux=mux;
    lmodes(ii).muy=muy;
    lmodes(ii).muz=muz;
    lmodes(ii).h=h;
    lmodes(ii).h_ = h_;
    
end
end

function out = matrixElements_q1_q2_q3_q4(n)

x1 = zeros(n,n);
x2 = zeros(n,n);
x3 = zeros(n,n);
x4 = zeros(n,n);

for i=0:n-1
    for j=0:n-1
        x1(i+1,j+1) = 1/sqrt(2)*(sqrt(j)*delta(i,j-1) + sqrt(j+1)*delta(i,j+1));
        x2(i+1,j+1) = (j+1/2) * delta(i,j) + ...
            1/2*sqrt((j+1)*(j+2)) * delta(i,j+2) ...
            +1/2*sqrt(j*(j-1)) * delta(i,j-2);
        x3(i+1,j+1) = 1/(2*sqrt(2))*(...
            sqrt(j*(j-1)*(j-2)) * delta(i,j-3)...
            + 3*j*sqrt(j)*delta(i,j-1)...
            + 3*(j+1)*sqrt(j+1) * delta(i,j+1)...
            + sqrt((j+1)*(j+2)*(j+3)) * delta(i,j+3));
        x4(i+1,j+1) = 1/4*(...
            sqrt(j*(j-1)*(j-2)*(j-3))* delta(i,j-4)...
            + 2*(2*j-1)*sqrt(j*(j-1)) * delta(i,j-2)...
            + (6*j^2+6*j+3) * delta(i,j)...
            + 2*(2*j+3)*sqrt((j+1)*(j+2)) * delta(i,j+2)...
            + sqrt((j+1)*(j+2)*(j+3)*(j+4))* delta(i,j+4));
    end
end

out.x1 = x1;
out.x2 = x2;
out.x3 = x3;
out.x4 = x4;
end

function out = delta(i,j)
if i==j
    out =1;
else
    out =0;
end
end