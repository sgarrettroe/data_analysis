function [V,E,ind] = analyzeEnergyLevels(lmodes,pmodes,varargin)

ind = 1:pmodes.NSTATES;
BW = realmax;
w_laser = 2000;
n_exciton_sig_figs = 1;
n_sparse_states = min(pmodes.NSTATES-2,200);

while length(varargin)>=2
    switch lower(varargin{1})
        case {'ind','state_list'}
            ind = varargin{2};
        case 'roptions'
            BW = varargin{2}.BW;
            w_laser = varargin{2}.w_laser;
        otherwise
            warning('unkonwn option %s\n',varargin{1})
    end
    
    varargin = varargin(3:end);
end

%upack the results
f=fieldnames(pmodes);
for ii=1:length(f)
    eval(strcat(f{ii},'=pmodes.',f{ii},';'))
end

% [original_energies,ordering] = sort(diag(H));
% original_gaps = original_energies - original_energies(1);
% original_vecs = IDENTITY(:,ordering);

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
    VV = speye(length(E),length(E)); %eigenvectors in eigenstate basis
else
    [original_vecs,original_energies] = eig(H,'vector');
    [original_energies,ordering] = sort(original_energies);
    original_vecs = original_vecs(:,ordering);
    % calculate the eigenvectors of the coupled system
    [V,E]=eig(H_,'vector');
    [E,ordering] = sort(E);
    V = V(:,ordering); %eigenvectors in input basis
    VV = eye(size(V)); %eigenvectors in eigenstate basis
end

new_energies = E;

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
for ii = iis;
    fprintf(1,'i = %3d\tE0 = %-8.1f E = %-8.1f gap0 = %-8.1f gap = %-8.1f \n',...
        ii,original_energies(ii),new_energies(ii),original_gaps(ii),new_gaps(ii));
end
fprintf(1,'\n');

%
iis = unique(vertcat(ind{:})');
for ii = iis;
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
    for ii = 1:nmodes
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
