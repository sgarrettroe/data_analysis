function [V,E] = analyzeEnergyLevels(lmodes,pmodes)

%upack the results
f=fieldnames(pmodes)
for ii=1:length(f)
    eval(strcat(f{ii},'=pmodes.',f{ii},';'))
end

% [original_energies,ordering] = sort(diag(H));
% original_gaps = original_energies - original_energies(1);
% original_vecs = IDENTITY(:,ordering);
[original_vecs,original_energies] = eig(H,'vector');
[original_energies,ordering] = sort(original_energies);
original_vecs = original_vecs(:,ordering);


original_gaps = original_energies - original_energies(1);

[V,E] = eig(H_,'vector');
new_energies = E;
new_gaps = new_energies - new_energies(1);
new_energies(2:end)-new_energies(1:end-1)

for ii = 1:min(6,length(V))
    fprintf(1,'i = %3d\tE0 = %-8.3f E = %-8.3f gap0 = %-8.3f gap = %-8.3f \n',...
        ii,original_energies(ii),new_energies(ii),original_gaps(ii),new_gaps(ii));
end
fprintf(1,'\n');

%
for ii = 1:min(6,length(V))
    fprintf(1,'i = %3d\tenergy gap = %8.3f\n', ii, new_gaps(ii));
    displayCoeffMatrix(V(:,ii),lmodes);
end

% transitions starting from ground state
fprintf('i\tj\to_gap\tn_gap\tu_orig\tu_mixed\n');
jjs = 1:min(6,length(V));
omu = zeros(1,length(jjs));
nmu = zeros(1,length(jjs));
for ii = 1
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

        fprintf(1,'%-8d%-8d%-8.3f%-8.3f%-8.3f%-8.3f\n',ii,jj,...
            original_gaps(jj),new_gaps(jj),omu(jj),nmu(jj));
        
    end
end
%fprintf(1,'\ndone\n\n');
