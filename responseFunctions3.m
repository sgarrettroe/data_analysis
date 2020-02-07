function out = responseFunctions3(pmodes,options)
global wavenumbersToInvPs k_B_SI h c_SI
k_B_cm_K = k_B_SI/h/c_SI/100;%k_B in cm-1/K 
thermal_cutoff = 0.01;
T = 0;
verbose = 0;

out = [];
order = options.order;

%canonical results
Rh.n = 2;
Rh.w = [1889.184 1947.621];
Rh.mu = [sqrt(0.908) 0 0 ; 0 sqrt(0.671) 0];
Rh.w_off = 1915;
Rh.w = Rh.w - Rh.w_off;
Rh.n2 = 3;
Rh.w2 = [3767.517 3883.911 3813.698];
Rh.w2 = Rh.w2 - 2*Rh.w_off;
Rh.mu2 = zeros(Rh.n2,3,Rh.n);
Rh.mu2(:,:,1) = [sqrt(2).*Rh.mu(1,:) ; 0,0,0 ; Rh.mu(2,:)];
Rh.mu2(:,:,2) = [ 0,0,0 ; sqrt(2).*Rh.mu(2,:) ; Rh.mu(1,:)];

% default values
flag_plot = false;
flag_print = true;
w0 = 0;
n_excitons = 2;
n_exciton_sig_figs = 1;

if isfield(options,'w0')
    if isempty(options.w0)
        w0=0;
    else
        w0 = options.w0;
    end
end
if isfield(options,'n_sparse_states')
    if isempty(options.n_sparse_states)
        %do nothing
    else
        n_sparse_states = options.n_sparse_states;
    end
end
if isfield(options,'flag_plot')
    flag_plot = options.flag_plot;
end
if flag_print
    fid = 1;
else
    if isunix
        fid = fopen('/dev/null');
    else
        warning('this windows "dev/null" is not tested yet');
        fid = fopen('nul');
    end        
end
if isfield(options,'T')
    if isempty(options.T)
        %do nothing (see default at top)
    else
        T = options.T;
    end
end
if isfield(options,'thermal_cutoff')
    if isempty(options.thermal_cutoff)
        %do nothing (see default at top)
    else
        thermal_cutoff = options.thermal_cutoff;
    end
end
if isfield(options,'verbose')
    if isempty(options.verbose)
        verbose=0;
    else
        verbose = options.w0;
    end
end

% simulation parameters
if verbose>=1,disp('initialize variables'),end
n_t = options.n_t;
n_zp = options.n_zp;
dt = options.dt; %time step
t2 = options.t2; %population time (ps)
pol = options.polarizations;
e_1 = pol{1};
e_2 = pol{2};
e_3 = pol{3};
e_4 = pol{4};
w_laser = options.w_laser;
BW = options.BW;

% set up time and response functions
J = zeros(1,n_t); 
J_accum = zeros(1,n_t);
R_r = zeros(n_t,n_t);
R_nr = zeros(n_t,n_t);
R_r_accum = zeros(n_t,n_t);
R_nr_accum = zeros(n_t,n_t);
t=0:dt:(n_t-1)*dt;
[T1,T3] = meshgrid(t,t);

% set up thermal weight
kT = k_B_cm_K*T;


%upack the results
if verbose>=1,disp('unpack pmodes'),end
f=fieldnames(pmodes);
for ii=1:length(f)
    eval(strcat(f{ii},'=pmodes.',f{ii},';'))
end

%
% set up operators
%
if issparse(H_)
    VV = speye(length(E),length(E)); %eigenvectors in eigenstate basis
else
    VV = eye(size(V)); %eigenvectors in eigenstate basis
end

if verbose>=1,disp('set up operators'),end
A0 = V'*A*V;
C0 = V'*C*V;

% rotate dipole operators to the eigenstate basis
MUX = V'*MUX*V;
MUY = V'*MUY*V;
MUZ = V'*MUZ*V;

% could add a loop over possible initial states here. The idea would be to
% look at the thermal density matrix elements relative to some cutoff. then
% loop through the response function calculation for each state with the
% appropriate thermal weight.

flag_finished_thermal_loop = false;

i_thermal = 0;
while ~flag_finished_thermal_loop
    i_thermal = i_thermal+1;
    if T == 0
        flag_finished_thermal_loop = true;
        thermal_weight = 1;
    else
        thermal_weight = exp(-E(i_thermal)/kT);
        if thermal_weight < thermal_cutoff
            flag_finished_thermal_loop = true;
            continue;
        end
    end
    if verbose>=1,fprintf(1,'loop over thermal states, i_thermal = %i\n',i_thermal);end
    if verbose>=1,fprintf(1,'    thermal weight = %f\n',thermal_weight);end
    
    if verbose>=1,disp('determine active states...'),end

% density matrix
PSIi = VV(:,i_thermal); %take first eigenstate for the time being
%rho = PSIi*PSIi'; % could do thermal density here!

%one way to go would be to define mui muj etc from inputs and then add
%invariants function (see thoughts below)
%tests:
%reproduce Rhcomplex 
%make sure in strong mixing dipoles are orthogonal
% check amplitudes of parallel and perp polarizations

% calc mus and omegas from inputs (ultimately want to refactor this)
% calculate the one and two exciton manifolds. Might need to be modified if
% thermal states are allowed. not sure. 

%find all one and two exciton states
[ind_1ex ind_2ex] =  findNExcitonStates(PSIi,C0,n_exciton_sig_figs);

%keep only the ones in the laser bandwidth
[ind_1ex ind_2ex] = filterExcitons(w_laser,BW,E,i_thermal,ind_1ex,ind_2ex);

% ind_1ex = find(abs(round(C0*PSIi,1))>0);
% ind_2ex = find(abs(round(C0*C0*PSIi,1))>0);
n = length(ind_1ex);
n2 = length(ind_2ex);

% energies -- subtract zero point energy
w = E(ind_1ex) - E(i_thermal);
w2 = E(ind_2ex) - E(i_thermal);

fprintf(fid,'one exciton state energies\n');
fprintf(fid,'%d\t%8.1f\n',[ind_1ex,w]');
fprintf(fid,'two exciton state energies\n');
fprintf(fid,'%d\t%8.1f\n',[ind_2ex,w2]');
fprintf(fid,'\n');

% subtract rotating frame frequency and convert to rad/ps
w = (w - w0)*2*pi*wavenumbersToInvPs;
w2 = (w2 - 2*w0)*2*pi*wavenumbersToInvPs;

if verbose>=1,disp('determine matrix elements...'),end

% calculate dipole matrix elements
mu = zeros(n,3);
mu2 = zeros(n2,3,n);
for ii = 1:n
    PSIf = VV(:,ind_1ex(ii));
    mu(ii,:) = [PSIf'*MUX*PSIi PSIf'*MUY*PSIi PSIf'*MUZ*PSIi];
    for jj =  1:n2
        PSIf2 = VV(:,ind_2ex(jj));
        mu2(jj,:,ii) = [PSIf2'*MUX*PSIf PSIf2'*MUY*PSIf PSIf2'*MUZ*PSIf];
    end
end

g = @(t) options.g(t,options.c2params);

if verbose>=1,disp('start linear spectroscopy...'),end
%linear spectroscopy
for j = 1:n
    [~,muj] = unit_vector(mu(j,:));
    J = J + muj^2.*exp(-1i*w(j).*t);
end
%add lineshape (same for all peaks for now)
J = J.*exp(-g(t));

if verbose>=1,disp('start third order spectroscopy...'),end
% first calculate all rephasing diagrams
fprintf(fid,'Rephasing transitions\n');
fprintf(fid,'i\tj\tk\tw_1\t(w_2)\tw_3\tu\n');
for j = 1:n
  for i  = 1:n
    
    [aa,mui] = unit_vector(mu(i,:));
    [bb,muj] = unit_vector(mu(j,:));
    dipole = mui^2*muj^2;
    angle = polarizationInvariant(e_1,e_2,e_3,e_4,...
        aa,bb,aa,bb);
    
    % rephasing diagram R1
    fprintf(fid,'%-8d%-8d%-8d%-8.1f%-8.1f%-8.1f%-8.3f\n',ind_1ex(j),ind_1ex(i),ind_1ex(j),w(j)./(2*pi*wavenumbersToInvPs)+w0,(w(j)-w(i))./(2*pi*wavenumbersToInvPs),w(i)./(2*pi*wavenumbersToInvPs)+w0,-dipole*angle);
    R_r = R_r - dipole*angle*exp(+ 1i*w(j).*T1 ...
				 - 1i*w(i).*T3 ...
				 + 1i*(w(j)-w(i))*t2);

             % rephasing diagram R2
    fprintf(fid,'%-8d%-8d%-8d%-8.1f%-8.1f%-8.1f%-8.3f\n',ind_1ex(j),ind_1ex(j),ind_1ex(i),w(j)./(2*pi*wavenumbersToInvPs)+w0,0,w(i)./(2*pi*wavenumbersToInvPs)+w0,-dipole*angle);
    R_r = R_r - dipole*angle*exp(+ 1i*w(j).*T1 ...
				 - 1i*w(i).*T3);
    
    for k = 1:n2      
      %molecular dipoles?
      [cc,muik_] = unit_vector(mu2(k,:,i));
      [dd,mujk_] = unit_vector(mu2(k,:,j));
      dipole = mui*muj*muik_*mujk_;
      angle = polarizationInvariant(e_1,e_2,e_3,e_4,...
         aa,bb,cc,dd);

     %rephasing diagram R3
       fprintf(fid,'%-8d%-8d%-8d%-8.1f%-8.1f%-8.1f%-8.3f\n',ind_1ex(j),ind_1ex(i),ind_2ex(k),w(j)./(2*pi*wavenumbersToInvPs)+w0,(w(j)-w(i))./(2*pi*wavenumbersToInvPs),(w2(k)-w(j))./(2*pi*wavenumbersToInvPs)+w0,dipole*angle);
       R_r = R_r + dipole*angle*exp(+ 1i*w(j).*T1 ...
				   - 1i*(w2(k)-w(j)).*T3 ...
				   + 1i*(w(j)-w(i)).*t2);
    end
  end
end
% add lineshape (same for all peaks for now)
R_r = exp(-g(T1)+g(t2)-g(T3)-g(T1+t2)-g(t2+T3)+g(T1+t2+T3)).*R_r;
fprintf(fid,'\n\n');

% now non-rephasing diagrams
fprintf(fid,'Non-rephasing transitions\n');
fprintf(fid,'i\tj\tk\tw_1\t(w_2)\tw_3\tu\n');
for j = 1:n
  for i  = 1:n      
    [aa,mui] = unit_vector(mu(i,:));
    [bb,muj] = unit_vector(mu(j,:));
    dipole = mui^2*muj^2;
    angle = polarizationInvariant(e_1,e_2,e_3,e_4,...
         aa,bb,aa,bb);
    
    % non-rephasing diagram R4
    fprintf(fid,'%-8d%-8d%-8d%-8.1f%-8.1f%-8.1f%-8.3f\n',ind_1ex(j),ind_1ex(i),ind_1ex(j),w(j)./(2*pi*wavenumbersToInvPs)+w0,(w(j)-w(i))./(2*pi*wavenumbersToInvPs),w(j)./(2*pi*wavenumbersToInvPs)+w0,-dipole*angle);
    R_nr = R_nr - dipole*angle*exp(- 1i*w(j).*T1 ...
				   - 1i*w(j).*T3 ... %?
				   - 1i*(w(j)-w(i))*t2);
    % non-rephasing diagram R5
    fprintf(fid,'%-8d%-8d%-8d%-8.1f%-8.1f%-8.1f%-8.3f\n',ind_1ex(j),ind_1ex(j),ind_1ex(i),w(j)./(2*pi*wavenumbersToInvPs)+w0,0,w(i)./(2*pi*wavenumbersToInvPs)+w0,-dipole*angle);
    R_nr = R_nr - dipole*angle*exp(- 1i*w(j).*T1 ...
				   - 1i*w(i).*T3);
    
    for k = 1:n2      
      %molecular dipoles
      [cc,muik_] = unit_vector(mu2(k,:,i));
      [dd,mujk_] = unit_vector(mu2(k,:,j));
      dipole = mui*muj*muik_*mujk_;
      angle = polarizationInvariant(e_1,e_2,e_3,e_4,...
         aa,bb,cc,dd);
     
      %non-rephasing diagram R6
      fprintf(fid,'%-8d%-8d%-8d%-8.1f%-8.1f%-8.1f%-8.3f\n',ind_1ex(j),ind_1ex(i),ind_2ex(k),w(j)./(2*pi*wavenumbersToInvPs)+w0,(w(j)-w(i))./(2*pi*wavenumbersToInvPs),(w2(k)-w(i))./(2*pi*wavenumbersToInvPs)+w0,dipole*angle);
      R_nr = R_nr + dipole*angle*exp(- 1i*w(j).*T1 ...
				     - 1i*(w2(k)-w(i)).*T3 ...
				     - 1i*(w(j)-w(i)).*t2);
    end
  end
end
% add lineshape (same for all peaks for now)
R_nr = exp(-g(T1)-g(t2)-g(T3)+g(T1+t2)+g(t2+T3)-g(T1+t2+T3)).*R_nr;

R_r_accum  = R_r_accum  + R_r *thermal_weight;
R_nr_accum = R_nr_accum + R_nr*thermal_weight;
J_accum = J_accum + J*thermal_weight;

fprintf(fid,'\n\n');
end

J = J_accum;
R_r = R_r_accum;
R_nr = R_nr_accum;

if verbose>=1,disp('calculate fourier transforms...'),end

% calculate 1D spectrum (freq domain)
J = real(fftshift(sgrsifft(J,n_zp)));

% divide first points (by row and column) by 2
R_r(:,1) = R_r(:,1)./2;
R_r(1,:) = R_r(1,:)./2;
R_nr(:,1) = R_nr(:,1)./2;
R_nr(1,:) = R_nr(1,:)./2;

if flag_plot
%what we have so far in the time domain
figure(1),clf
subplot(1,2,1)
contourf(real(R_r'),10); 
axis equal tight
subplot(1,2,2)
contourf(real(R_nr'),10);
axis equal tight
end

% do the fft
R_r = ifft2(R_r,n_zp,n_zp); %given the frequency definitions used
                            %above, use the ifft to get the
                            %frequencies right (Mathematica has the
                            %opposite definition of the fft by default)
R_nr = ifft2(R_nr,n_zp,n_zp);

%this is the frequency not the energy of the transitions
freq  = fftFreqAxis(t,'time_units','ps','zeropad',n_zp);
freq = freq+w0;

if flag_plot
%now frequency domain
figure(2),clf
subplot(1,2,1)
contourf(freq,freq,fftshift(real(R_r')),20); %pump-probe axis convention
%contourf(fftshift(real(R_r)),20; % the (omega_1, omega_3) axis convention
axis equal tight
subplot(1,2,2)
contourf(freq,freq,fftshift(real(R_nr')),20)
%contourf(fftshift(real(R_nr)),20)
axis equal tight
end

% flip R_r (being careful to keep zero frequency as the first time
% point), add the response functions, take the real part, and
% finally reorganize so that the 0 frequency is in the center
R = fftshift(real(fliplr(circshift(R_r,[0 -1]))+R_nr));
    

if flag_plot
figure(3),clf
n_contours = 40;
MAX = max(abs(R(:)));
level_list = linspace(-MAX,MAX,n_contours+2);
dl = level_list(2)-level_list(1);
cmin =level_list(1)-dl/2;
cmax =level_list(end);

%contourf(freq,freq,R',level_list) %use R' to display the pump-probe axis convention
contourf(freq,freq,R,level_list) %use R to display the (omega_1, omega_3) axis convention
caxis([cmin cmax]);			
axis equal tight
end

if verbose>=1,disp('package output...'),end

%package output
out.w1 = freq;
out.w3 = freq;
out.J = J;
out.R_r = R_r;
out.R_nr = R_nr;
out.R = R;
out.E = E;
out.V = V;
out.ind_1ex = ind_1ex;
out.ind_2ex = ind_2ex;
out.energy_gap1 = w./(2*pi*wavenumbersToInvPs)+w0;
out.energy_gap2 = w2./(2*pi*wavenumbersToInvPs)+2*w0;
out.mu1 = mu;
out.mu2 = mu2;

end

function [out,n] = unit_vector(in)
n  = norm(in,2);
if n>1e-6
out = in./n;
else
    n=0;
    out = 0*in;
end
out = out(:); %convert to a column matrix
end

function n_sparse_states = estimateNSparseStates(pmodes,options)
original_energies = diag(pmodes.H);
original_energies = original_energies-original_energies(1);
max_e = 2*(options.w_laser+options.BW);
n_sparse_states = sum(original_energies<=max_e);
end