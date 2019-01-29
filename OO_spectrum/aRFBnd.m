classdef (Abstract) aRFBnd < aRF
    
    
    properties
        lb; %lower bounds
        ub; %upper bounds
        useParallel=false;
    end
    
    methods
        function obj = aRFBnd(options)
            % assign properties from the input options struct
            props = properties(obj);
            for ii = 1:length(props)
                if isfield(options,props{ii})
                    obj.(props{ii}) = options.(props{ii});
                end
            end
            
            % get parameters initialized
            obj = obj.makeParamStruct;
            obj.freeParamNames = obj.freeFitParamNames;
            obj.p0 = obj.freeFitParamInitialValues;
            obj.lb = obj.freeFitParamLowerBounds;
            obj.ub = obj.freeFitParamUpperBounds;
            
            %move to super-class?
            obj = obj.setupTimeAxes; %has to come before freq
            obj = obj.setupFreqAxes;
            obj = obj.setupSimMatrix;
            %obj = obj.calcSpectrum;
            
            %for subclasses?
            %obj = obj.addNoise;
            %obj = obj.bandwidthFiltering;
            %obj = ...
        end
        
        function [out] = freeFitParamInitialValues(obj)
            props = properties(obj);
            out = [];
            for ii = 1:length(props)
                if isa(obj.(props{ii}),'fitParam')
                    if obj.(props{ii}).isFree
                        out = [out,obj.(props{ii}).value];
                    end
                end
            end
        end
        function [out] = freeFitParamLowerBounds(obj)
            props = properties(obj);
            out = [];
            for ii = 1:length(props)
                if isa(obj.(props{ii}),'fitParam')
                    if obj.(props{ii}).isFree
                        out = [out,obj.(props{ii}).lb];
                    end
                end
            end
        end
        function [out] = freeFitParamUpperBounds(obj)
            props = properties(obj);
            out = [];
            for ii = 1:length(props)
                if isa(obj.(props{ii}),'fitParam')
                    if obj.(props{ii}).isFree
                        out = [out,obj.(props{ii}).ub];
                    end
                end
            end
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
            [pfit,err] = fmincon(@(p)obj.err_fun(p),obj.p0,[],[],[],[],obj.lb,obj.ub,[],opt);
            
            toc
            
            obj.pfit = pfit;
            obj.err = err;
            
            for ii = 1:length(pfit)
                fprintf(1,'%20s\t%12f\n',obj.freeParamNames{ii},pfit(ii));
            end
        end
    end
end
