classdef multiMeasurement
    % class that can serve as a holder for multiple measurements (like
    % polarizations, and adds the chi2 of each in the error calculation
    %
    % the first spectrum sets the initial guess and the fitting tolerances
    properties
        s; %holds the aRFs
        p0;
        pfit; %output fit
        err; % final error
        useParallel;
        nboot;
        pboot;
        tolfun;
        tolx;
        maxfun;
        lb;
        ub;
        freeParamNames;
    end
    methods
        function obj = multiMeasurement(s_in)
            obj.s = s_in;
            obj.p0 = s_in(1).p0;
            obj.nboot = s_in(1).nboot;
            obj.useParallel = s_in(1).useParallel;
            obj.tolfun = s_in(1).tolfun;
            obj.tolx = s_in(1).tolx;
            obj.maxfun = s_in(1).maxfun;
            obj.lb = s_in(1).lb;
            obj.ub = s_in(1).ub;
            obj.freeParamNames = s_in(1).freeParamNames;
        end
        
        function obj = globalFit(obj)%dataMatrix,w1,w3,p_array,gfstruct,lb,ub)
            
            % set fit parameters (I believe based on version of matlab)
            if exist('optimoptions','file'),
                opt = optimoptions('fmincon','Algorithm','active-set','Display','iter','TolFun',obj.tolfun,'TolX',obj.tolx,'MaxFunEvals',obj.maxfun);
            elseif exist('optimset','file'),
                %    opt = optimset('Algorithm','active-set','Display','iter','TolFun',tolfun,'TolX',tolx,'MaxFunEvals',maxfun);
                opt = optimset('Algorithm','active-set','Display','iter','TolFun',obj.tolfun,'TolX',obj.tolx,'MaxFunEvals',obj.maxfun);
            else
                warning('Could not set parameters! Look for the right options functon for your installation.');
            end
            
            opt = optimoptions(opt,'UseParallel',obj.useParallel);

            tic
            
            %fmincon has required parameters of error function and initial guess.
            %Documentation has a bunch of additional parameters, most of which we don't
            %understand, but the syntax for not using them is to leave them as blanks.
            [obj.pfit,obj.err] = fmincon(@(p)obj.err_fun(p),obj.p0,[],[],[],[],obj.lb,obj.ub,[],opt);
            
            toc
            
            for ii = 1:length(obj.pfit)
                fprintf(1,'%20s\t%12f\n',obj.freeParamNames{ii},obj.pfit(ii));
            end
        end
        
        function obj = globalFitBootstrap(obj)%dataMatrix,w1,w3,p_array,gfstruct,lb,ub)
            
            % set fit parameters (I believe based on version of matlab)
            if exist('optimoptions','file'),
                opt = optimoptions('fmincon','Algorithm','active-set','Display','iter','TolFun',obj.tolfun,'TolX',obj.tolx,'MaxFunEvals',obj.maxfun);
            elseif exist('optimset','file'),
                %    opt = optimset('Algorithm','active-set','Display','iter','TolFun',tolfun,'TolX',tolx,'MaxFunEvals',maxfun);
                opt = optimset('Algorithm','active-set','Display','iter','TolFun',obj.tolfun,'TolX',obj.tolx,'MaxFunEvals',obj.maxfun);
            else
                warning('Could not set parameters! Look for the right options functon for your installation.');
            end
            
            opt = optimoptions(opt,'UseParallel',obj.useParallel);
            
            obj.pboot = zeros(obj.nboot,length(obj.p0));
            
            npoints = zeros(1,length(obj.s));
            for ii = 1:length(obj.s)
                npoints(ii) = numel(obj.s(ii).dataMatrix);
            end
            ind = cell(1,length(obj.s));
            
            % go!
            bootStart = tic;
            for ii = 1:obj.nboot;
                
                for jj = 1:length(obj.s)
                    ind{jj} = randi(npoints(jj),1,npoints(jj));
                end
                
                %[pboot(ii,:),err] = fmincon(@(p)obj.err_fun(p),obj.p0,[],[],[],[],obj.lb,obj.ub,[],opt);
                [obj.pboot(ii,:)] = fmincon(@(p,ind)obj.err_fun_boot(p,ind),obj.p0,[],[],[],[],obj.lb,obj.ub,[],opt,ind);
            end
            bootEnd = toc(bootStart);
            
            fprintf(1,'Elapsed time: %i seconds\n',bootEnd);
            % the 2 sigma limits should be 95% intervals
            for ii = 1:size(obj.pboot,2)
                fprintf(1,'%20s\t%12f\tpm\t%12.3f\n',...
                    obj.s(1).freeParamNames{ii},mean(obj.pboot(:,ii)),2*std(obj.pboot(:,ii)));
            end
            
        end
        
        function chi2 = err_fun(obj,p)
            chi2 = 0;
            for ii = 1:length(obj.s)
                chi2 = chi2+obj.s(ii).err_fun(p);
            end
        end
        function chi2 = err_fun_boot(obj,p,ind)
            chi2 = 0;
            for ii = 1:length(obj.s)
                chi2 = chi2+obj.s(ii).err_fun_boot(p,ind{ii});
            end
            
        end
        
    end
end
