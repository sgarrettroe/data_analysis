function out = responseFunctions(pmodes,lmodes,options)
% calculate response functions
% TODO
% 1. dipole tensor products with isotropic averaging
% 2. g(t) lineshape functions
%X3. speed
%    ? sparse matrices (around 25 states sparse starts to be faster...)
%    X preeval or functions for A(t) C(t) 
%X4. recast PSIi*PSIi' as rho 
% 4b. recast rho = exp(-beta H_)/Trace(exp(-beta H_) to add
%     thermal distribution (ie for CO2 bend modes)
% 5. add center frequency properly --> interaction representation?
%X6. move params to input structure
%X7. the cross-peaks still look funny. what is going on? --> in Rhcomplex.m
%    too, so I guess it is right...
%X8. add energy cutoff or something to limit number of states in the calc?
% 9. "vectorize" the for loops

out = [];
order = options.order;

% default values
flag_plot = false;
flag_sparse = false;
flag_energy_cutoff = false;
w0 = 0;

if isfield(options,'flag_sparse')
    flag_sparse = logical(options.flag_sparse);
end
if isfield(options,'energy_cutoff')
    if isempty(options.energy_cutoff)
        flag_energy_cutoff = false;
        energy_cutoff = options.energy_cutoff;
    else
        flag_energy_cutoff = true;
        energy_cutoff = options.energy_cutoff;
    end
    
    if isfield(options,'flag_energy_cutoff')
        %in case you want to explicitly turn off the energy cutoff but
        %specify a value for the energy_cutoff anyway
        flag_energy_cutoff = options.flag_energy_cutoff;
    end
end
if isfield(options,'w0')
    if isempty(options.w0)
        w0=0;
    else
        w0 = options.w0;
    end
end
if isfield(options,'flag_plot')
    flag_plot = options.flag_plot;
end

if w0~=0,
    error('the energy offset does not work yet!!! set 0 or empty to proceed.'),
end

% simulation parameters
n_t = options.n_t;
n_zp = options.n_zp;
dt = options.dt; %time step
t2 = options.t2; %population time (ps)
T2 = options.T2; %dephasing time (ps)


%upack the results
f=fieldnames(pmodes);
for ii=1:length(f)
    eval(strcat(f{ii},'=pmodes.',f{ii},';'))
end

[V,E]=eig(H_,'vector');
[E,ordering] = sort(E);
V = V(:,ordering);

%
% set up operators
%

% reduce size of the active space
if flag_energy_cutoff
    % this is assuming things are sorted by energy
    
    % creation and annihil ops in sorted order
    A0 = V'*A*V;
    C0 = V'*C*V;
    
    % find the states below the cutoff
    ind = ((E-E(1))<=energy_cutoff);
    
    % keep lower energies only
    E = E(ind);
   
    E = E*2*pi*wavenumbersToInvPs; %convert cm-1 to rad/ps

    % extract those parts of the creation, annhil, and eigenvector ops
    A0 = A0(ind,ind);
    C0 = C0(ind,ind);
    V = V(ind,ind);
else
    %if not using the energy cutoff
    
    E = E*2*pi*2.9979e10/1e12; %convert cm-1 to rad/ps

    A0 = V'*A*V;
    C0 = V'*C*V;
end

% time shifted creation C and annihilation A operators
A = @(t) diag(exp(1i*E*t))*A0*diag(exp(-1i*E*t));
C = @(t) diag(exp(1i*E*t))*C0*diag(exp(-1i*E*t));

% density matrix
PSIi = V(:,1); %take lowest energy state
rho = PSIi*PSIi'; % could do thermal density here!

if flag_sparse
    rho = sparse(rho);
    A0 = sparse(A0);
    C0 = sparse(C0);
    A = @(t) sparse(1:NSTATES,1:NSTATES,exp(1i*E*t))...
        *A0*sparse(1:NSTATES,1:NSTATES,exp(-1i*E*t));
    C = @(t) sparse(1:NSTATES,1:NSTATES,exp(1i*E*t))...
        *C0*sparse(1:NSTATES,1:NSTATES,exp(-1i*E*t));
end    

% 
% Parameters of the simulation

t=0:dt:(n_t-1)*dt;
[T1,T3] = meshgrid(t,t);

w = fftFreqAxis(t,'time_units','ps','zeropad',n_zp);
w = w+w0;
time = t;

% linear response
if order>=1
J = zeros(1,n_t);

for ii = 1:n_t
    t1 = t(ii);
%    J(ii) = PSIi'*expm(1i*H_*(t1))*A*expm(-1i*H_*(t1))*C*PSIi;
    J(ii) = trace(A(t1)*C0*rho);
end

J = J.*exp(-t/T2);
J_w = fftshift(sgrsifft(J,n_zp));
    
if flag_plot
    figure(1),clf
    plot(t,real(J),t,imag(J))
    figure(2),clf
    plot(w,real(J_w),w,imag(J_w))
end

out.t = time;
out.w = w;
out.J = J;
out.J_w = J_w;
end

% nonlinear response with creation annihilation ops
if order>=3
R_r = zeros(n_t,n_t);
R_nr = zeros(n_t,n_t);

%     mui = sqrt(mu(i,:)*mu(i,:)');
%     muj = sqrt(mu(j,:)*mu(j,:)');
%     cos1 = mu(i,:)*mu(j,:)'/(mui*muj);
%     angle = (1 + 2*cos1^2 )/15;
%     dipole = mui^2*muj^2;
%     % rephasing diagram R1
%     R_r = R_r - dipole*angle*exp(+ 1i*w(j)*c.*T1 ...
% 				 - 1i*w(i)*c.*T3 ...
% 				 + 1i*(w(j)-w(i))*c*t2 ...
% 				 - (T1+T3)/T2);
%     % rephasing diagram R2
%     R_r = R_r - dipole*angle*exp(+ 1i*w(j)*c.*T1 ...
% 				 - 1i*w(i)*c.*T3 ...
% 				 - (T1+T3)/T2);
tic
for ii = 1:n_t
    for jj = 1:n_t
        t1 = t(ii);
        t3 = t(jj);
       % GSB
       %<~00
       %->10
       %<-00
       %->10
%         R_nr(jj,ii) = R_nr(jj,ii) - trace(expm(1i*H_*(t1+t2+t3))*A*expm(-1i*H_*(t1+t2+t3))...
%             *expm(1i*H_*(t1+t2))*C*expm(-1i*H_*(t1+t2))...
%             *expm(1i*H_*(t1))*A*expm(-1i*H_*(t1))...
%             *C*(PSIi*PSIi'));      
        R_nr(jj,ii) = R_nr(jj,ii) ...
            - trace(A(t1+t2+t3)*C(t1+t2)*A(t1)*C0*rho);      
       % SE
       %<~00
       %  10->
       %  11<-
       %->10
%        R_nr(jj,ii) = R_nr(jj,ii) - trace(expm(1i*H_*(t1+t2+t3))*A*expm(-1i*H_*(t1+t2+t3))...
%            *C*(PSIi*PSIi')...
%            *expm(1i*H_*(t1))*A*expm(-1i*H_*(t1))...
%            *expm(1i*H_*(t1+t2))*C*expm(-1i*H_*(t1+t2))...
%            );
       R_nr(jj,ii) = R_nr(jj,ii) ...
           - trace(A(t1+t2+t3)*C0*rho*A(t1)*C(t1+t2));
        % GSB
        %<~00   t1+t2+t3
        %->10   t1+t2
        %  00-> t1
        %  01<- t0
%         R_r(jj,ii) = R_r(jj,ii) - trace(expm(1i*H_*(t1+t2+t3))*A*expm(-1i*H_*(t1+t2+t3))...
%             *expm(1i*H_*(t1+t2))*C*expm(-1i*H_*(t1+t2))...
%             *(PSIi*PSIi')...
%             *A*expm(1i*H_*(t1))*C*expm(-1i*H_*(t1)));   
        R_r(jj,ii) = R_r(jj,ii) ...
            - trace(A(t1+t2+t3)*C(t1+t2)*rho*A0*C(t1));   
        % SE
        %<~00   t1+t2+t3
        %  10-> t1+t2
        %->11   t1
        %  01<- t0
%         R_r(jj,ii) = R_r(jj,ii) - trace(expm(1i*H_*(t1+t2+t3))*A*expm(-1i*H_*(t1+t2+t3))...
%             *expm(1i*H_*(t1))*C*expm(-1i*H_*(t1))...
%             *(PSIi*PSIi')...
%             *A...
%             *expm(1i*H_*(t1+t2))*C*expm(-1i*H_*(t1+t2))...
%             );
        R_r(jj,ii) = R_r(jj,ii) ...
            - trace(A(t1+t2+t3)*C(t1)*rho*A0*C(t1+t2));
        % ESA
        %<~00   t1+t2+t3
        %->21   t1+t2
        %  11<- t1
        %->10   t0
%         R_nr(jj,ii) = R_nr(jj,ii) + trace(expm(1i*H_*(t1+t2+t3))*A*expm(-1i*H_*(t1+t2+t3))...
%             *expm(1i*H_*(t1+t2))*C*expm(-1i*H_*(t1+t2))...
%             *C*(PSIi*PSIi')...
%             *expm(1i*H_*(t1))*A*expm(-1i*H_*(t1)));      
        R_nr(jj,ii) = R_nr(jj,ii) ...
            + trace(A(t1+t2+t3)*C(t1+t2)*C0*rho*A(t1));      
        % ESA
        %<~11   t1+t2+t3
        %->21   t1+t2
        %->11   t1
        %  01<- t0
%         R_r(jj,ii) = R_r(jj,ii) + trace(expm(1i*H_*(t1+t2+t3))*A*expm(-1i*H_*(t1+t2+t3))...
%             *expm(1i*H_*(t1+t2))*C*expm(-1i*H_*(t1+t2))...
%             *expm(1i*H_*(t1))*C*expm(-1i*H_*(t1))...
%             *(PSIi*PSIi')...
%             *A);   
        R_r(jj,ii) = R_r(jj,ii) ...
            + trace(A(t1+t2+t3)*C(t1+t2)*C(t1)*rho*A0);   
        
    end
end
toc

R_r = R_r.*exp(-(T1+T3)/T2);
R_nr = R_nr.*exp(-(T1+T3)/T2);

if flag_plot
    %what we have so far in the time domain
    figure(11),clf
    subplot(1,2,1)
    contourf(real(R_r'),10);
    axis equal tight
    subplot(1,2,2)
    contourf(real(R_nr'),10);
    axis equal tight
end

% ffts
R_r  = sgrsifft2(R_r, n_zp);
R_nr = sgrsifft2(R_nr,n_zp);

%
%now frequency domain

if flag_plot
    figure(3),clf
    subplot(1,2,1)
    title('rephasing')
    %contourf(fftshift(real(R_r'A)),20); %pump-probe axis convention
    contourf(w,w,fftshift(real(R_r)),20); % the (omega_1, omega_3) axis convention
    axis equal tight
    subplot(1,2,2)
    title('non-rephasing')
    %contourf(w,w,fftshift(real(R_nr')),20)
    contourf(w,w,fftshift(real(R_nr)),20)
    axis equal tight
end

%
% flip R_r (being careful to keep zero frequency as the first time
% point), add the response functions, take the real part, and
% finally reorganize so that the 0 frequency is in the center
R = fftshift(real(fliplr(circshift(R_r,[0 -1]))+R_nr));
    
if flag_plot
    figure(4),clf
    n_contours = 40;
    MAX = max(abs(R(:)));
    level_list = linspace(-MAX,MAX,n_contours+2);
    dl = level_list(2)-level_list(1);
    cmin =level_list(1)-dl/2;
    cmax =level_list(end);
    
    %contourf(w,w,R,level_list) %use R' to display the pump-probe axis convention
    contourf(w,w,R,level_list) %use R to display the (omega_1, omega_3) axis convention
    caxis([cmin cmax]);
    axis equal tight
end

out.R_r = R_r;
out.R_nr = R_nr;
out.R = R;
end

out.E = E;
out.V = V;
