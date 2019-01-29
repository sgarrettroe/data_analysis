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
    end
    methods
        function obj = multiMeasurement(s_in)
            obj.s = s_in;
            obj.p0 = s_in(1).p0;
        end
        
        function obj = globalFit(obj)%dataMatrix,w1,w3,p_array,gfstruct,lb,ub)
            
            % set fit parameters (I believe based on version of matlab)
            if exist('optimoptions','file'),
                opt = optimoptions('fmincon','Algorithm','active-set','Display','iter','TolFun',obj.s(1).tolfun,'TolX',obj.s(1).tolx,'MaxFunEvals',obj.s(1).maxfun);
            elseif exist('optimset','file'),
                %    opt = optimset('Algorithm','active-set','Display','iter','TolFun',tolfun,'TolX',tolx,'MaxFunEvals',maxfun);
                opt = optimset('Algorithm','active-set','Display','iter','TolFun',obj.s(1).tolfun,'TolX',obj.s(1).tolx,'MaxFunEvals',obj.s(1).maxfun);
            else
                warning('Could not set parameters! Look for the right options functon for your installation.');
            end
            
            
            tic
            
            %fmincon has required parameters of error function and initial guess.
            %Documentation has a bunch of additional parameters, most of which we don't
            %understand, but the syntax for not using them is to leave them as blanks.
            [obj.pfit,obj.err] = fmincon(@(p)obj.err_fun(p),obj.s(1).p0,[],[],[],[],obj.s(1).lb,obj.s(1).ub,[],opt);
            
            toc
            
            for ii = 1:length(obj.pfit)
                fprintf(1,'%20s\t%12f\n',obj.s(1).freeParamNames{ii},obj.pfit(ii));
            end
        end
        
        function chi2 = err_fun(obj,p)
            chi2 = 0;
            for ii = 1:length(obj.s)
                chi2 = chi2+obj.s(ii).err_fun(p);
            end
        end
        
    end
end
