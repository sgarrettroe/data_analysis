function result = correlationTime(cfit_obj,varargin)
%function [tau_c,stde_tau_c] = correlationTime(cfit_obj,varargin)
%
% v1 Tom Brinzer 2017-10-04
% v2 SGR 2023-05-10

% parse inputs
default_t_initial = 0;
default_t_final = Inf;

p = inputParser;
validCFitObj = @(x) isa(x,'cfit') | ...
    (iscell(x) && all(cellfun(@(xx)isa(xx,'cfit'),x)));
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validScalarNonNegNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
addRequired(p,'cfit_obj',validCFitObj);
addParameter(p,'t_initial',default_t_initial,validScalarNonNegNum);
addParameter(p,'t_final',default_t_final,validScalarPosNum);
parse(p,cfit_obj,varargin{:});

cfit_obj  = p.Results.cfit_obj;
t_initial = p.Results.t_initial;
t_final   = p.Results.t_final;

if ~iscell(cfit_obj)
    cfit_obj = {cfit_obj};
end

result = correlationTimeBackend(cfit_obj,t_initial,t_final);

end

function result = correlationTimeBackend(cfit_cell,t_initial,t_final)
n_models = length(cfit_cell);
result(n_models) = struct('tau_c',[],'stde_tau_c',[]);

for ii = 1:n_models
    cfit_obj = cfit_cell{ii};
    
    % find the locations of our amplitudes and times in the fit for
    % later (this 'should' be consistent, but I'd rather check, since
    % that'd be a disgusting thing to mess up).
    names = coeffnames(cfit_obj);
    
    % Complain if there is a parameter that looks like an offset (b or c)
    bTest = regexp(names,'[bc]\d*');
    ind3 = ~cellfun(@isempty,bTest);
    if any([bTest{:}]) && t_final==Inf
        warning('Fit function appears to contain a constant offset and upper limit is Inf! Results may not be accurate with an infinite upper limit. Consider restricting the upper limit, t_final, to the length of the data.')
    end
    
    % times
    tTest = regexp(names,'t[12345]');
    ind1 = ~cellfun(@isempty,tTest);
    
    % amplitudes
    aTest = regexp(names,'a[12345]');
    ind2 = ~cellfun(@isempty,aTest);
    ind = ind1 + ind2 + ind3;
    if sum(~ind) > 0
        n = join(names);
        error(['Variable naming format is incorrect. Expected coefficients '...
            'with names a1 t1 etc; got %s'],n{:})
    end
    
    
    % calculate the correlation time
    tau_c = integrate(cfit_obj,t_final,t_initial);
    
    % calculate the standard deviation of the correlation time:
    % make a convenient matrix of coefficient values and standard
    % errors
    values = coeffvalues(cfit_obj);
    CI = confint(cfit_obj);
    stde = (CI(2,:)-CI(1,:))./(2*1.96);
    coeffMat = [values' stde'];
    
    % and then decomposing it into its time and amplitude components
    % I suppose this code will break if amplitudes and times are not in
    % the same order ... something to keep in mind for the future
    tMat = coeffMat(ind1,:);
    aMat = coeffMat(ind2,:);
    
    % Propagation of error from the best fit coefficients
    stde_tau_c = sqrt(sum(tMat(:,1).^2.*aMat(:,2).^2 + ...
        aMat(:,1).^2.*tMat(:,2).^2));
    
    % package results
    result(ii).tau_c = tau_c;
    result(ii).stde_tau_c = stde_tau_c;
    
end
end