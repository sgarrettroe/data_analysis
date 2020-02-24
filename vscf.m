function scfmodes = vscf(lmodes,tod,fod,varargin)
% try to do some kind of vibrational SCF calculation to make a basis that
% includes anharmonicity in a mean field sense. Hopefully this will reduce
% the size of the space needed for the VCI calculations (diagonalizing the
% vibrational hamiltonian in the product mode basis

energy_convergence_criterion=0.01; %cm-1
flag_print = 1;
flag_pause = 0;
if flag_print
    fid=1;
else
    if isunix
        fid=fopen('/dev/null','w');
    elseif ispc
        fid = fopen('nul','w');
    end
end

scfmodes = lmodes;
n_modes = length(lmodes);
for ii = 1:n_modes;
    scfmodes(ii).U = lmodes(ii).identity;
end
scfmodes_old = scfmodes;


% active_modes = 1:length(lmodes);
% if nargin == 4
%     active_modes = varargin{1};
% end

%start scf loop
i_scf = 0;
flag_done = false;
while ~flag_done
    i_scf = i_scf + 1;
    %start loop over modes
    for i_mode = 1:n_modes;
        F = zeros(size(lmodes(i_mode).identity));
        
        order = 3;
        n_tod = size(tod,1);
        for ii = 1:n_tod
            inds = tod(ii,1:order);
            val = tod(ii,order+1);
            
            % break loop if no mode is our active mode
            if ~ismember(i_mode,inds)
                continue
            end
            
            F_i = calcNextTerm(i_mode,inds,val,scfmodes_old);
            
            F = F + F_i;
        end
        
        
        order = 4;
        n_fod = size(fod,1);
        for ii = 1:n_fod
            inds = fod(ii,1:order);
            val = fod(ii,order+1);
            % break loop if no mode is our active mode
            if ~ismember(i_mode,inds)
                continue
            end
            
            F_i = calcNextTerm(i_mode,inds,val,scfmodes_old);
            
            F = F + F_i;
            
            %             n_perms = nperms(inds);
            %
            %             ind_active = (inds==i_mode);
            %
            %             OP = cell(1,4);
            %             for jj = 1:length(OP)
            %
            %                 if ind_active(jj)
            %                     OP{jj} = lmodes(inds(jj)).q;
            %                 else
            %                     v = scfmodes_old(inds(jj)).U(:,1);
            %                     OP{jj} = v'*lmodes(inds(jj)).q*v;
            %                 end
            %             end
            %
            %
            %             F = F + 1/24*n_perms*val*OP{1}*OP{2}*OP{3}*OP{4};
            
        end
        
        % add single mode (harmonic) hamiltonian
        F = F + scfmodes_old(i_mode).h;
        S = scfmodes_old(i_mode).U' * scfmodes_old(i_mode).U;
        [V,D] = eig(F,S); %Vs are orthogonal but not normalized...
        [D,ordering] = sort(diag(D));
        V = V(:,ordering);
        norms = sqrt(diag(V'*V));
        scfmodes(i_mode).U = V./repmat(norms',[length(norms) 1]);
        scfmodes(i_mode).E = D;
        
        if flag_print
            n_plot_points = length(D);
            x_axis = 1:n_plot_points;
            y_axis1 = diag(scfmodes(i_mode).h);
            y_axis2 = D;
            y_axis3 = y_axis1(:)-y_axis2(:);
            figure(1)
            subplot(n_modes,1,i_mode),plot(x_axis,y_axis1,'-o',x_axis,y_axis2,x_axis,y_axis3);
            title(['difference between harmonic and scf energies for mode ',num2str(i_mode),' step ',num2str(i_scf)]);
            figure(2)
            if i_scf==1&&i_mode==1,clf,hold on,end
            plot(x_axis,y_axis3)
        end
        
    end %i_modes
    
    delta_E = zeros(1,n_modes);
    if i_scf>1
        for i_modes = 1:n_modes
            delta_E(i_modes) = scfmodes(i_modes).E(1) - scfmodes_old(i_modes).E(1);
        end
        if all(abs(delta_E)<=energy_convergence_criterion)
            flag_done = true;
        end
    end
    
    if flag_print && i_scf==1
        fprintf(fid,'i_scf \tdEs\n');
    end
    fprintf(fid,'%i\t',i_scf);
    for i_mode = 1:n_modes
        fprintf(fid,'%f\t',delta_E(i_mode));
    end
    fprintf(fid,'\n');
    
    
    if ~flag_done&&flag_pause,pause,end
    
    %update the coefficients of the basis vectors for the next pass
    scfmodes_old = scfmodes;
    
end %i_scf

%when the calc has converged, rotate all ops to the new basis
ops = {'a','c','q','p','mux','muy','muz','h','h_'};
n_ops = length(ops);
for i_modes = 1:n_modes
    
    U = scfmodes(i_modes).U;
    
    for i_ops = 1:n_ops
        
        O = scfmodes(i_modes).(ops{i_ops});
        scfmodes(i_modes).(ops{i_ops}) = U'*O*U;
        
    end
    
end
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

function  F_i = calcNextTerm(i_mode,inds,val,scfmodes_old)
[n_syms,n_unique,unique_syms,num_repeats] = countSymbols(inds);
%n_perms = nperms(n_syms,n_repeats);
%prefactor = 1/factorial(n_syms);
prefactor = 1/prod(factorial(num_repeats));
[unique_syms,n] = sortSymbols(n_unique,unique_syms,num_repeats,i_mode);

OP = scfmodes_old(i_mode).q^n(1);
for i_loop = 2:n_unique
    
    tmp = unique_syms(i_loop);
    v = scfmodes_old(tmp).U(:,1);
    aux = v'*scfmodes_old(tmp).q^n(i_loop)*v;
    
    OP = OP*aux;
end
%F_i = prefactor*n_perms*val*OP;
F_i = prefactor*val*OP;
end