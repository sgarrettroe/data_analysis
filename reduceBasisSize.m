function modes = reduceBasisSize(modes,inds_cell)

n_modes = length(modes);
if n_modes~=length(inds_cell)
    error('size mismatch between n_modes=%i and input cell=%i',n_modes,length(inds_cell))
end

for i_mode = 1:n_modes
    inds = inds_cell{i_mode};
    
    n_keep = length(inds);
    n_states = size(modes(i_mode).identity,1);
    
    % projector
    P = zeros(n_states,n_keep);
    for ii = 1:length(inds)
        P(inds(ii),ii) = 1;
    end
    
    modes(i_mode).nstates = length(inds);
    if isfield(modes,'E')
        modes(i_mode).E = P'*modes(i_mode).E;
    end
    
    ops = {'identity','a','c','q','p','mux','muy','muz','h','h_','U','q_pow_1','q_pow_2','q_pow_3','q_pow_4'};
    n_ops = length(ops);
    for ii = 1:n_ops
        
        if isfield(modes,ops{ii})
            modes(i_mode).(ops{ii}) = P'*modes(i_mode).(ops{ii})*P;
        end
        
    end
    
    %renormalize the basis states U
    if isfield(modes,'U')
        norms = sqrt(diag(modes(i_mode).U'*modes(i_mode).U));
        modes(i_mode).U = modes(i_mode).U./repmat(norms',[length(norms) 1]);
    end
end
