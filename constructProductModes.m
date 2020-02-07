function pmodes = constructProductModes(lmodes)
% constructProductModes takes an array of local modes as an input and
% builds from that a tensor product of operators and basis functions to
% span the joint space.

% some general things
pmodes.NSTATES = prod([lmodes.nstates]);
pmodes.IDENTITY = eye(pmodes.NSTATES);

% number of modes we have to deal with
n_lmodes = length(lmodes);

%build product modes out of ops defined for each local mode
ops = {'a','c','q','mux','muy','muz','h'};
%ops = fieldnames(lmodes);
for ii = 1:length(ops)
    % set up big cell array to hold all the operators that we will be
    % tensor multiplying
    big_cell = cell(n_lmodes,n_lmodes);
    for jj = 1:n_lmodes
        %copy the identities everywhere
        for kk = 1:n_lmodes
            big_cell{kk,jj} = lmodes(jj).identity; %put identities in first
        end
        big_cell{jj,jj} = lmodes(jj).(ops{ii}); %now add the operator itself
    end
    
    for jj = 1:n_lmodes
        %resulting operator
        OP = cumkron(big_cell(jj,:));
        
        %make a new field with the name of the input field and the mode number
        %all made uppercase
        fname = upper(strcat(ops{ii},num2str(jj)));
        pmodes.(fname) = OP;
    end
    
end

%make total operators out of the modes
ops = {'MUX','MUY','MUZ','H','A','C'};
for ii = 1:length(ops)
    temp = zeros(size(pmodes.IDENTITY));
    for jj = 1:n_lmodes
        %field name
        fname = strcat(ops{ii},num2str(jj));
        %accumulate
        temp = temp + pmodes.(fname);
    end
    pmodes.(upper(ops{ii})) = temp;
end

%make empty H_ for Hamiltonian with couplings (to be added later)
pmodes.H_ = [];

end
