function pmodes = reduceProductBasisSize(pmodes,energy_cutoff)

E = full(diag(pmodes.H));
inds = find(E<=energy_cutoff);

    n_keep = length(inds);
    n_states = size(pmodes.H,1);
    
    % projector
    P = zeros(n_states,n_keep);
    for ii = 1:length(inds)
        P(inds(ii),ii) = 1;
    end
    
    pmodes.NSTATES = length(inds);
    if isfield(pmodes,'E')
        pmodes.E = P'*pmodes.E;
    end
    
%    ops = {'identity','a','c','q','p','mux','muy','muz','h','h_','U','q_pow_1','q_pow_2','q_pow_3','q_pow_4'};
ops = {'H','H_','C','A','MUX','MUY','MUZ'};
    n_ops = length(ops);
    for ii = 1:n_ops
        
        if isfield(pmodes,ops{ii})
            pmodes.(ops{ii}) = full(P'*pmodes.(ops{ii})*P);
        end
        
    end
    
    % throw away the rest of the ops that we don't need
    rmops = fieldnames(pmodes);
    for i=1:length(rmops)
        if any(strcmpi(rmops{i},ops))||strcmpi(rmops{i},'NSTATES')
            %don't erase what we need
            continue
        else
            %throw away the rest
            pmodes.(rmops{i}) = [];
        end
    end
    
    %renormalize the basis states U
    if isfield(pmodes,'U')
        norms = sqrt(diag(pmodes(i_mode).U'*pmodes(i_mode).U));
        pmodes(i_mode).U = pmodes(i_mode).U./repmat(norms',[length(norms) 1]);
    end
end
