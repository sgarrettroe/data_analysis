function [pboot,gfstruct] = gfBootstrap_w(dataMatrix,WtMatrix,w1,w3,p_array,gfstruct,lb,ub,nboot)

% Rather than changing our data matrix intensity, we now scale our response
% function calculation to it.
scale = gfScale(dataMatrix);
maxScale = abs(min(min(min(dataMatrix))));
% err_fun = @(p,gfstruct) sum(sum(sum(((dataMatrix-analyticalResponseFunctionsFun(p,w1,w3,gfstruct)).*WtMatrix).^2)));
err_fun = @(p,gfstruct) sum(sum(sum(((dataMatrix(gfstruct.bootstrap)-analyticalResponseFunctionsFun(p,w1,w3,gfstruct).*scale(gfstruct.bootstrap)).^2).*WtMatrix(gfstruct.bootstrap))))./maxScale;

% 
% err_fun_b = @(p,gfstruct) sum(sum(sum(((dataMatrix(gfstruct.bootstrap) - analyticalResponseFunctionsFun(p,w1,w3,gfstruct)).*WtMatrix(gfstruct.bootstrap)).^2.)));


pboot = zeros(nboot,length(p_array));
npoints = numel(dataMatrix);

tolfun = npoints*1e-15;
tolx = 1e-5;
maxfun = 1e3;

% set fit parameters
if exist('optimoptions')==2,
    opt = optimoptions('fmincon','Algorithm','active-set','Display','iter','TolFun',tolfun,'TolX',tolx,'MaxFunEvals',maxfun);
elseif exist('optimset')==2,
    %    opt = optimset('Algorithm','active-set','Display','iter','TolFun',tolfun,'TolX',tolx,'MaxFunEvals',maxfun);
    opt = optimset('Algorithm','active-set','Display','iter','TolFun',tolfun,'TolX',tolx,'MaxFunEvals',maxfun);
else
    warning('Could not set parameters! Look for the right options functon for your installation.');
end

% go!
bootStart = tic;
parfor ii = 1:nboot;
    tempglobalfit = struct(gfstruct);
    tempglobalfit.bootstrap = randi(npoints,1,npoints);
    
    [pboot(ii,:)] = fmincon(err_fun,p_array,[],[],[],[],lb,ub,[],opt,tempglobalfit);
end
bootEnd = toc(bootStart);

fprintf(1,'Elapsed time: %i seconds\n',bootEnd);
% the 2 sigma limits should be 95% intervals
for ii = 1:size(pboot,2)
    fprintf(1,'%20s\t%12f\tpm\t%12.3f\n',...
        gfstruct.pnames{ii},mean(pboot(:,ii)),2*std(pboot(:,ii)));
end

% 
% %show the distribution of taus
% figure,hist(pboot(:,2))
% title('Histogram of tau from bootstrap analysis');
% xlabel('tau / ps')
% ylabel('counts')