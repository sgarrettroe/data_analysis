function [pfit,err,gfstruct] = globalFit_w(dataMatrix,WtMatrix,w1,w3,p_array,gfstruct,lb,ub)


% Rather than changing our data matrix intensity, we now scale our response
% function calculation to it.
scale = gfScale(dataMatrix);
maxScale = abs(min(min(min(dataMatrix))));
% err_fun = @(p,gfstruct) sum(sum(sum(((dataMatrix-analyticalResponseFunctionsFun(p,w1,w3,gfstruct)).*WtMatrix).^2)));
err_fun = @(p,gfstruct) sum(sum(sum(((dataMatrix-analyticalResponseFunctionsFun(p,w1,w3,gfstruct).*scale).^2).*WtMatrix)))./maxScale;


%tolfun relates to the squared residual in the error
%tolx relates to the displacement in your parameters
tolfun = 1e-17;
tolx = 1e-10;
maxfun = 1e3; %maximum number of function evaluations


% set fit parameters (I believe based on version of matlab)
if exist('optimoptions')==2,
    opt = optimoptions('fmincon','Algorithm','active-set','Display','iter','TolFun',tolfun,'TolX',tolx,'MaxFunEvals',maxfun);
elseif exist('optimset')==2,
%    opt = optimset('Algorithm','active-set','Display','iter','TolFun',tolfun,'TolX',tolx,'MaxFunEvals',maxfun);
    opt = optimset('Algorithm','active-set','Display','iter','TolFun',tolfun,'TolX',tolx,'MaxFunEvals',maxfun);
else
    warning('Could not set parameters! Look for the right options functon for your installation.');
end


tic

%fmincon has required parameters of error function and initial guess.
%Documentation has a bunch of additional parameters, most of which we don't
%understand, but the syntax for not using them is to leave them as blanks.

[pfit,err] = fmincon(err_fun,p_array,[],[],[],[],lb,ub,[],opt,gfstruct)
toc

gfstruct.pfit = pfit;
gfstruct.fitMethod = 'weighted fit (individual points)';

for ii = 1:length(pfit)
    fprintf(1,'%20s\t%12f\n',gfstruct.pnames{ii},pfit(ii))
end