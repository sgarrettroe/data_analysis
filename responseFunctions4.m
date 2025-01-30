function out = responseFunctions4(V, E, pmodes, varargin)
%responseFunctions4 Calculate 2D spectra from VCI and VPT2

global wavenumbersToInvPs k_B_SI h c_SI
k_B_cm_K = k_B_SI/h/c_SI/100;%k_B in cm-1/K

% default values
default_BW = 500;
default_c2form = '1fast';
default_c2params = struct('T2', 3);
default_dt = 0.025;  % time step
default_flag_interstate_coherences = true;
default_flag_orientational_response = false;
default_flag_plot = false;
default_flag_print = true;
default_n_exciton_tol = 0.0001;
default_n_t = 128;
default_n_zp = 256;  % 2 * nt
default_T = 0;
default_t2 = 0;  % population time (ps)
default_tau_orient = inf;
default_thermal_cutoff = 0.01;
default_verbose = 0;
default_w0 = 0;
default_w_laser = 2000;
default_polarizations = {[1; 0; 0], [1; 0; 0], [1; 0; 0], [1; 0; 0] };

out = [];

p = inputParser;
addRequired(p, 'V')
addRequired(p, 'E')
addRequired(p, 'pmodes')
addParameter(p, 'BW', default_BW, @(x) isnumeric(x) & isscalar(x));
addParameter(p, 'c2form', default_c2form, @(x) ischar(x));
addParameter(p, 'c2params', default_c2params, @(x) isstruct(x));
addParameter(p, 'dt', default_dt, @(x) isnumeric(x) & isscalar(x));
addParameter(p, 'flag_interstate_coherences', ...
    default_flag_interstate_coherences, @(x) islogical(x) & isscalar(x));
addParameter(p, 'flag_orientational_response', ...
    default_flag_orientational_response, ...
    @(x) islogical(x) & isscalar(x));
addParameter(p, 'flag_plot', default_flag_plot, ...
    @(x) islogical(x) & isscalar(x));
addParameter(p, 'flag_print', default_flag_print, ...
    @(x) islogical(x) & isscalar(x));
addParameter(p, 'n_exciton_tol', default_n_exciton_tol, ...
    @(x) isnumeric(x) & isscalar(x));
addParameter(p, 'n_t', default_n_t, @(x) isnumeric(x) & isscalar(x));
addParameter(p, 'n_zp', default_n_zp, @(x) isnumeric(x) & isscalar(x));
addParameter(p, 'order', @(x) x==3);
addParameter(p, 'T', default_T, ...
    @(x) isnumeric(x) & isscalar(x));
addParameter(p, 'polarizations', default_polarizations, @(x) iscell(x) & numel(x)==4);
addParameter(p, 't2', default_t2, @(x) isnumeric(x) & isscalar(x));
addParameter(p, 'tau_orient', default_tau_orient, ...
    @(x) isnumeric(x) & isscalar(x));
addParameter(p, 'thermal_cutoff', default_thermal_cutoff, ...
    @(x) isnumeric(x) & isscalar(x));
addParameter(p, 'verbose', default_verbose, ...
        @(x) isnumeric(x) & isscalar(x));
addParameter(p, 'w0', default_w0, ...
    @(x) isnumeric(x) & isscalar(x));
addParameter(p, 'w_laser', default_w_laser, ...
    @(x) isnumeric(x) & isscalar(x));

parse(p, V, E, pmodes, varargin{:});

BW = p.Results.BW;
c2form = p.Results.c2form;
c2params = p.Results.c2params;
dt = p.Results.dt;
flag_interstate_coherences = p.Results.flag_interstate_coherences;
flag_orientational_response = p.Results.flag_orientational_response;
flag_plot = p.Results.flag_plot;
flag_print = p.Results.flag_print;
n_exciton_tol = p.Results.n_exciton_tol;
n_t = p.Results.n_t;
n_zp = p.Results.n_zp;
T  = p.Results.T;
t2  = p.Results.t2;
tau_orient = p.Results.tau_orient;
thermal_cutoff  = p.Results.thermal_cutoff;
verbose  = p.Results.verbose;
w0 = p.Results.w0;
w_laser = p.Results.w_laser;
polarizations = p.Results.polarizations;
pol = polarizations;

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

% simulation parameters
if verbose>=1,disp('initialize variables'),end

% set up lineshape function
lineshape = chooseLineshapeFunction(c2form, c2params);
g = @(t) lineshape.g(t, c2params);

% set up polarization
e_1 = pol{1};
e_2 = pol{2};
e_3 = pol{3};
e_4 = pol{4};
if all(e_1==e_2) && all(e_1==e_3) && all(e_1==e_4)
    polarization_experiment_name = 'parallel';
elseif all(e_1==e_2) && e_1'*e_3==0 && all(e_3==e_4)
    polarization_experiment_name = 'perpendicular';
elseif (all(e_1==e_3) && e_1'*e_2==0 && all(e_2==e_4) )|...
        (all(e_1==e_4) && e_1'*e_2==0 && all(e_2==e_3))
    polarization_experiment_name = 'perpendicular';
else
    warning('unknown polarization scheme, orientational response turned OFF')
    flag_orientational_response = false;
end

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
MUX = pmodes.MUX;
MUY = pmodes.MUY;
MUZ = pmodes.MUZ;
H_ = pmodes.H_;

%
% Energies and eigenvectors are now inputs
%

% VV is eigenvectors in eigenstate basis
if issparse(H_)
    VV = speye(length(E),length(E));
else
    VV = eye(size(V));
end

%
% set up operators
%
if verbose>=1,disp('set up operators'),end

% rotate dipole operators to the eigenstate basis
MUX = V'*MUX*V;
MUY = V'*MUY*V;
MUZ = V'*MUZ*V;

flag_finished_thermal_loop = false;
i_thermal = 0;
while ~flag_finished_thermal_loop
    i_thermal = i_thermal + 1;
    if T == 0
        flag_finished_thermal_loop = true;
        thermal_weight = 1;
    else
        thermal_weight = exp( - ( E(i_thermal) - E(1) ) / kT );
        if thermal_weight < thermal_cutoff
            if i_thermal ~= 1
                flag_finished_thermal_loop = true;
                continue;
            end
        end
    end
    if verbose>=1,fprintf(1,'loop over thermal states, i_thermal = %i\n',i_thermal);end
    if verbose>=1,fprintf(1,'    thermal weight = %f\n',thermal_weight);end

    if verbose>=1,disp('determine active states...'),end

    % density matrix
    PSIi = VV(:,i_thermal); %take first eigenstate for the time being

    lower_limit_energy = w_laser - BW/2;
    upper_limit_energy = w_laser + BW/2;

    % calculate dipole matrix elements
    n = length(E);
    mu = zeros(n,3);
    mu2 = zeros(n,3,n);
    for ii = i_thermal+1:n-1
        PSIf = VV(:,ii);
        Ef = E(ii)-E(i_thermal);

        %filter out energies outside our window
        if (Ef < lower_limit_energy)||(Ef > upper_limit_energy)
            continue;
        end

        mu(ii,:) = [PSIf'*MUX*PSIi, PSIf'*MUY*PSIi, PSIf'*MUZ*PSIi];
        for jj =  (ii+1):n
            PSIf2 = VV(:,jj);

            Ef2 = E(jj)-Ef-E(i_thermal);

            %filter out energies outside our window
            if (Ef2 < lower_limit_energy)||(Ef2 > upper_limit_energy)
                continue;
            end

            mu2(jj,:,ii) = [PSIf2'*MUX*PSIf, PSIf2'*MUY*PSIf, PSIf2'*MUZ*PSIf];
        end
    end


    ind_1ex = find(sum(mu.^2,2) > n_exciton_tol);
    ind_2ex = [];
    for ii = 1:length(ind_1ex)
        ind_2ex = [ind_2ex; ...
            find(sum(mu2(:,:,ind_1ex(ii)).^2,2) > n_exciton_tol)];
    end
    ind_2ex = unique(ind_2ex);

    %reduce dipole matrix elements to only the needed size
    mu = mu(ind_1ex,:);
    mu2 = mu2(ind_2ex,:,ind_1ex);

    n = length(ind_1ex);
    n2 = length(ind_2ex);

    % energies -- subtract zero point energy
    w = E(ind_1ex) - E(i_thermal);
    w2 = E(ind_2ex) - E(i_thermal);

    fprintf(fid, 'one exciton state energies\n');
    fprintf(fid, '%d\t%8.1f\n', [ind_1ex,w]');
    fprintf(fid, 'two exciton state energies\n');
    fprintf(fid, '%d\t%8.1f\n', [ind_2ex,w2]');
    fprintf(fid, '\n');

    % subtract rotating frame frequency and convert to rad/ps
    w = (w - w0) * 2 * pi * wavenumbersToInvPs;
    w2 = (w2 - 2 * w0) * 2 * pi * wavenumbersToInvPs;

    if verbose>=1,disp('determine matrix elements...'),end

    if verbose>=1,disp('start linear spectroscopy...'),end
    %linear spectroscopy
    for j = 1:n
        [~,muj] = unit_vector(mu(j,:));
        J = J + muj^2 .* exp(-1i * w(j) .* t);
    end
    %add lineshape (same for all peaks for now)
    J = J .* exp(-g(t));

    if verbose>=1, disp('start third order spectroscopy...'), end
    % first calculate all rephasing diagrams
    fprintf(fid, 'Rephasing transitions\n');
    fprintf(fid, 'i\tj\tk\tw_1\t(w_2)\tw_3\tu\n');
    for j = 1:n
        if flag_interstate_coherences
            loop_indices = 1:n;
        else
            loop_indices = j;
        end
        for i  = loop_indices

            [aa, mui] = unit_vector(mu(i,:));
            [bb, muj] = unit_vector(mu(j,:));
            dipole = mui^2 * muj^2;
            angle = polarizationInvariant(e_1, e_2, e_3, e_4,...
                aa, bb, aa, bb);

            % rephasing diagram R1
            fprintf(fid,'%-8d%-8d%-8d%-8.1f%-8.1f%-8.1f%-8.3f\n',ind_1ex(j),ind_1ex(i),ind_1ex(j),w(j)./(2*pi*wavenumbersToInvPs)+w0,(w(j)-w(i))./(2*pi*wavenumbersToInvPs),w(i)./(2*pi*wavenumbersToInvPs)+w0,-dipole*angle);
            R_r = R_r - dipole * angle * exp(+ 1i * w(j) .* T1 ...
                - 1i * w(i) .* T3 ...
                + 1i * (w(j) - w(i)) * t2);

            % rephasing diagram R2
            fprintf(fid,'%-8d%-8d%-8d%-8.1f%-8.1f%-8.1f%-8.3f\n',ind_1ex(j),ind_1ex(j),ind_1ex(i),w(j)./(2*pi*wavenumbersToInvPs)+w0,0,w(i)./(2*pi*wavenumbersToInvPs)+w0,-dipole*angle);
            R_r = R_r - dipole * angle * exp(+ 1i * w(j) .* T1 ...
                - 1i * w(i) .* T3);

            for k = 1:n2
                %molecular dipoles?
                [cc, muik_] = unit_vector(mu2(k,:,i));
                [dd, mujk_] = unit_vector(mu2(k,:,j));
                dipole = mui * muj * muik_ * mujk_;
                angle = polarizationInvariant(e_1, e_2, e_3, e_4,...
                    aa, bb, cc, dd);

                %rephasing diagram R3
                fprintf(fid,'%-8d%-8d%-8d%-8.1f%-8.1f%-8.1f%-8.3f\n',ind_1ex(j),ind_1ex(i),ind_2ex(k),w(j)./(2*pi*wavenumbersToInvPs)+w0,(w(j)-w(i))./(2*pi*wavenumbersToInvPs),(w2(k)-w(j))./(2*pi*wavenumbersToInvPs)+w0,dipole*angle);
                R_r = R_r + dipole * angle * exp(+ 1i * w(j) .* T1 ...
 				   - 1i * (w2(k)-w(j)) .* T3 ...
                   + 1i * (w(j)-w(i)) .* t2);
            end
        end
    end
    % add lineshape (same for all peaks for now)
    R_r = exp(-g(T1) + g(t2) - g(T3) - g(T1 + t2) - g(t2 + T3) + g(T1 + t2 + T3)) .* R_r;

    % add orientational relaxation (same for all peaks now)
    if flag_orientational_response
        [p, q, r] = orientationalResponse(tau_orient, 3, T1, t2, T3);
        if strcmpi(polarization_experiment_name, 'parallel')
            R_r = R_r .* p;
        elseif strcmpi(polarization_experiment_name, 'perpendicular')
            R_r = R_r .* q;
        elseif strcmpi(polarization_experiment_name, 'crossed')
            R_r = R_r .* r;
        end
    end
    fprintf(fid, '\n\n');

    % now non-rephasing diagrams
    fprintf(fid, 'Non-rephasing transitions\n');
    fprintf(fid, 'i\tj\tk\tw_1\t(w_2)\tw_3\tu\n');
    for j = 1:n
        if flag_interstate_coherences
            loop_indices = 1:n;
        else
            loop_indices = j;
        end
        for i  = loop_indices
            [aa, mui] = unit_vector(mu(i,:));
            [bb, muj] = unit_vector(mu(j,:));
            dipole = mui^2 * muj^2;
            angle = polarizationInvariant(e_1, e_2, e_3, e_4, ...
                aa, bb, aa, bb);

            % non-rephasing diagram R4
            fprintf(fid,'%-8d%-8d%-8d%-8.1f%-8.1f%-8.1f%-8.3f\n',ind_1ex(j),ind_1ex(i),ind_1ex(j),w(j)./(2*pi*wavenumbersToInvPs)+w0,(w(j)-w(i))./(2*pi*wavenumbersToInvPs),w(j)./(2*pi*wavenumbersToInvPs)+w0,-dipole*angle);
            R_nr = R_nr - dipole * angle * exp(- 1i * w(j) .* T1 ...
                - 1i * w(j) .* T3 ...
                - 1i * (w(j)-w(i)) * t2);
            % non-rephasing diagram R5
            fprintf(fid,'%-8d%-8d%-8d%-8.1f%-8.1f%-8.1f%-8.3f\n',ind_1ex(j),ind_1ex(j),ind_1ex(i),w(j)./(2*pi*wavenumbersToInvPs)+w0,0,w(i)./(2*pi*wavenumbersToInvPs)+w0,-dipole*angle);
            R_nr = R_nr - dipole * angle * exp(-1i * w(j) .* T1 ...
                - 1i * w(i) .* T3);

            for k = 1:n2
                %molecular dipoles
                [cc, muik_] = unit_vector(mu2(k,:,i));
                [dd, mujk_] = unit_vector(mu2(k,:,j));
                dipole = mui * muj * muik_ * mujk_;
                angle = polarizationInvariant(e_1, e_2, e_3, e_4,...
                    aa, bb, cc, dd);

                %non-rephasing diagram R6
                fprintf(fid,'%-8d%-8d%-8d%-8.1f%-8.1f%-8.1f%-8.3f\n',ind_1ex(j),ind_1ex(i),ind_2ex(k),w(j)./(2*pi*wavenumbersToInvPs)+w0,(w(j)-w(i))./(2*pi*wavenumbersToInvPs),(w2(k)-w(i))./(2*pi*wavenumbersToInvPs)+w0,dipole*angle);
                R_nr = R_nr + dipole * angle * exp(-1i * w(j) .* T1 ...
                    - 1i * (w2(k)-w(i)) .* T3 ...
                    - 1i * (w(j)-w(i)) .* t2);
            end
        end
    end
    % add lineshape (same for all peaks for now)
    R_nr = exp(-g(T1) - g(t2) - g(T3) + g(T1 + t2) + g(t2 + T3) - g(T1 + t2 + T3)) .* R_nr;

    % add orientational relaxation (same for all peaks now)
    if flag_orientational_response
        [p, q, r] = orientationalResponse(tau_orient, 3, T1, t2, T3);
        if strcmpi(polarization_experiment_name,'parallel')
            R_nr = R_nr .* p;
        elseif strcmpi(polarization_experiment_name, 'perpendicular')
            R_nr = R_nr .* q;
        elseif strcmpi(polarization_experiment_name, 'crossed')
            R_nr = R_nr .* r;
        end
    end
    fprintf(fid,'\n\n');

    R_r_accum  = R_r_accum  + R_r * thermal_weight;
    R_nr_accum = R_nr_accum + R_nr * thermal_weight;
    J_accum = J_accum + J * thermal_weight;

    fprintf(fid,'\n\n');
end

J = J_accum;
R_r = R_r_accum;
R_nr = R_nr_accum;

if verbose>=1, disp('calculate fourier transforms...'), end

% calculate 1D spectrum (freq domain)
J = real(fftshift(sgrsifft(J,n_zp)));

% divide first points (by row and column) by 2
R_r(:,1) = R_r(:,1) ./ 2;
R_r(1,:) = R_r(1,:) ./ 2;
R_nr(:,1) = R_nr(:,1) ./ 2;
R_nr(1,:) = R_nr(1,:) ./ 2;

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
    clim([cmin cmax]);
    axis equal tight
end

if verbose>=1, disp('package output...'), end

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