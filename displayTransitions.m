function displayTransitions(response,lmodes,pmodes)

ind_1ex = response.ind_1ex;
ind_2ex = response.ind_2ex;

n = length(ind_1ex);
n2 = length(ind_2ex);

% energy differences w/o perturbation -- subtract zero point energy
E = diag(pmodes.H);
[E,ordering] = sort(E);
VV = pmodes.IDENTITY; %eigenvectors in original basis
VV = VV(:,ordering); %eigenvectors in input basis

original_gaps = E(ind_1ex) - E(1);
original_gaps2 = E(ind_2ex) - E(1);

% calculate dipole matrix elements
PSIi = VV(:,1); %take first eigenstate for the time being

MUX = pmodes.MUX;
MUY = pmodes.MUY;
MUZ = pmodes.MUZ;

original_mu = zeros(n,3);
original_mu2 = zeros(n2,3,n);
for ii = 1:n
    PSIf = VV(:,ind_1ex(ii));
    original_mu(ii,:) = [PSIf'*MUX*PSIi PSIf'*MUY*PSIi PSIf'*MUZ*PSIi];
    for jj =  1:n2
        PSIf2 = VV(:,ind_2ex(jj));
        original_mu2(jj,:,ii) = [PSIf2'*MUX*PSIf PSIf2'*MUY*PSIf PSIf2'*MUZ*PSIf];
    end
end

V = response.V;
new_gaps = response.energy_gap1;
new_gaps2 = response.energy_gap2;
new_mu = response.mu1;
new_mu2 = response.mu2;

% display energy levels of involved transitions
for ii = 1:n
    fprintf(1,'i = %3d\to_gap = %-8.1f n_gap = %-8.1f \n',...
        ii,original_gaps(ii),new_gaps(ii));
end
fprintf(1,'\n');
for ii = 1:n2
    fprintf(1,'j = %3d\to_gap = %-8.1f n_gap = %-8.1f \n',...
        ind_2ex(ii),original_gaps2(ii),new_gaps2(ii));
end
fprintf(1,'\n');

disp('ONE EXCITON STATES');
% display coefficient matrix for each involved state
for ii = 1:n
    fprintf(1,'j = %3d\to_gap = %8.3f\tn_gap = %8.3f\n', ind_1ex(ii), original_gaps(ii),new_gaps(ii));
    displayCoeffMatrix(V(:,ind_1ex(ii)),lmodes);
end
disp('TWO EXCITON STATES');
for ii = 1:n2
    fprintf(1,'i = %3d\to_gap = %8.3f\tn_gap = %8.3f\n', ii, original_gaps2(ii),new_gaps2(ii));
    displayCoeffMatrix(V(:,ind_2ex(ii)),lmodes);
end

% transitions starting from ground state
for ii = 1:n
    fprintf('i\tj\to_gap\tn_gap\tu_orig\tu_mixed\n');
    for jj = 1:n2
        omu = (original_mu(ii,:)*original_mu(ii,:)')...
            *(original_mu2(jj,:,ii)*original_mu2(jj,:,ii)');
        nmu = (new_mu(ii,:)*new_mu(ii,:)')...
            *(new_mu2(jj,:,ii)*new_mu2(jj,:,ii)');
        
        fprintf(1,'%-8d%-8d%-8.1f%-8.1f%-8.2f%-8.2f\n',ind_1ex(ii),ind_2ex(jj),...
            original_gaps2(jj),new_gaps2(jj),omu,nmu);
        
    end
    fprintf(1,'\n\n');
end
%fprintf(1,'\ndone\n\n');
