%% try bootstrapping the constrained optimization with weights
pboot = zeros(nboot,length(p0));
npoints = numel(temp.data);

tolfun = npoints*1e-15;
tolx = 1e-5;
maxfun = 1e3;
% set fit parameters
if exist('optimoptions')==2,
    opt = optimoptions('fmincon','Algorithm','active-set','Display','iter','TolFun',tolfun,'TolX',tolx,'MaxFunEvals',maxfun);
elseif exist('optimset')==2,
    opt = optimset('Algorithm','active-set','Display','iter','TolFun',tolfun,'TolX',tolx,'MaxFunEvals',maxfun);
else
    warning('Could not set parameters! Look for the right options functon for your installation.');
end

% go!
tic
parfor ii = 1:nboot;
    disp(ii)
    tempglobalfit = struct(globalfit);
    tempglobalfit.bootstrap = randi(npoints,1,npoints);
    
    [pboot(ii,:)] = fmincon(temp.err_fun_b,p0,[],[],[],[],lb,ub,[],opt,tempglobalfit);
    %originally, this fmincon was working on 'globalfit,' however, the
    %variable is changing in the loops, in parallel with itself ... so we
    %defined a new dummy variable that then is allowed to change
end
toc

globalfit.bootstrap = pboot;
%%
% the 2 sigma limits should be 95% intervals
for ii = 1:size(pboot,2)
    fprintf(1,'%20s\t%12f\tpm\t%12.3f\n',...
        globalfit.pnames{ii},mean(globalfit.bootstrap(:,ii)),2*std(globalfit.bootstrap(:,ii)));
end

%%

beep
flag_save = input('Save output? (1 = "YES"; 0 = "NO"): ');
if flag_save
    save(sprintf('%s-%02i_summary.mat',datestring{1},experiment{1}),...
        'globalfit','-append')
end
flag_save = 0;