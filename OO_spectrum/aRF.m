classdef (Abstract) aRF
    
    
    properties
        %
        %  Things that can't be fit parameters
        %
        
        %time axis properties
        dt;
        n_t;
        n_zp;
        n_under = 0;
        t;
        
        
        % time axes
        T1;
        t2;
        T3;
        
        % frequency axes
        w;
        W1;
        W3;
        dw;
        flag_rotating_frame = true; %make hidden?
        
        
        % experimental
        w1_in;
        w3_in;
        %polarization;%deprecated
        t2_array;
        dataMatrix;
        
        % response functions for each feynmanDiagram
        diagrams;
        n_diagrams;
        
        % fit parameters
        %tolfun relates to the squared residual in the error
        %tolx relates to the displacement in your parameters
        tolfun = 1e-6;%1e-17?
        tolx = 1e-6;%1e-10?
        maxfun = 1e3; %maximum number of function evaluations
        
        %container for fit parameters, these are the ones actually used for fitting
        paramStruct; %structure of all parameters (fixed and free)
        freeParamNames; %names of only free params (cell array)
        p0; %starting point (vector)
        weightMatrix; %optional for _w functions
        
        %
        % Things that might be fit parameters
        %
        
        %center frequency
        w_01_cm;
        %phase
        phase_deg;
        
        %lineshapeFunction
        damping;
        
        %
        %  Output
        %
        simMatrix;
        spec; %a spectrum just for testing
        pfit; %ending point (vector)
        err; %resulting error
    end
    
    methods (Abstract)
        makeResponseFunctions(obj,p);%deprecated, should clean up!!!
        calcSpectrum(obj,p);
    end
    methods
        function obj = aRF(options)
            if nargin>0
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
                
                %move to super-class?
                obj = obj.setupTimeAxes; %has to come before freq
                obj = obj.setupFreqAxes;
                obj = obj.setupSimMatrix;
                %obj = obj.calcSpectrum;
                
            end
        end
        
        function obj = updateFreeFitParams(obj,p)
            for ii = 1:length(obj.freeParamNames)
                obj.paramStruct.(obj.freeParamNames{ii}) = p(ii);
            end
        end
        
        function obj = makeParamStruct(obj)
            n = obj.allFitParamNames;
            v = num2cell(obj.allFitParamInitialValues);
            obj.paramStruct = cell2struct(v,n,2);
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
        function [out] = freeFitParamValues(obj)
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
        function [out] = freeFitParamNames(obj)
            props = properties(obj);
            out = {};
            for ii = 1:length(props)
                if isa(obj.(props{ii}),'fitParam')
                    for jj = 1:length(obj.(props{ii}))
                        if obj.(props{ii})(jj).isFree
                            n = obj.(props{ii})(jj).name;
                            out = [out,n]; %square brackets make it so cells don't get nested
                        end
                        %what about out = [out n(obj.(props{ii}).isFree]; if
                        %isFree is made to be an array
                    end
                end
            end
        end
        function [out] = allFitParamInitialValues(obj)
            props = properties(obj);
            out = [];
            for ii = 1:length(props)
                if isa(obj.(props{ii}),'fitParam')
                    
                    out = [out,obj.(props{ii}).value];
                    
                end
            end
        end
        function [out] = allFitParamNames(obj)
            props = properties(obj);
            out = {};
            for ii = 1:length(props)
                if isa(obj.(props{ii}),'fitParam')
                    n = [obj.(props{ii}).name];
                    out = [out,n]; %the square brackets make it so the cells don't get nested
                end
            end
        end
        
        function obj = calcDiagramsTime(obj)
            for jj=1:obj.n_diagrams
                obj.diagrams(jj) = obj.diagrams(jj).calcResponseTime(obj.T1,obj.t2,obj.T3);
            end
        end
        
        function obj = calcDiagramsFreq(obj,p)
            for jj=1:obj.n_diagrams
                obj.diagrams(jj) = obj.diagrams(jj).timeToFreq(p);
            end
        end
        
        function obj = addDiagrams(obj)
            obj.spec = zeros(obj.n_zp,obj.n_zp);
            for ii = 1:obj.n_diagrams
                obj.spec = obj.spec + obj.diagrams(ii).R;
            end
        end
        
        function obj = resample(obj,i_t2)
            %obj.spec = obj.spec./abs(min(obj.spec(:))); %normalize to the 01 band (negative)
            obj.simMatrix(:,:,i_t2) = interp2(obj.W1,obj.W3,obj.spec,obj.w1_in,obj.w3_in','*linear');
        end
        
        function obj = setupFreqAxes(obj)
            
            obj.w = fftFreqAxis(obj.t,...
                'time_units','ps',...
                'freq','wavenumbers',...
                'shift','on',...
                'zeropad',obj.n_zp,...
                'undersampling',obj.n_under);
            if obj.flag_rotating_frame
                obj.w = obj.w + obj.paramStruct.w_01_cm;
            end
            obj.dw = obj.w(2)-obj.w(1);
            [obj.W1,obj.W3] = meshgrid(obj.w,obj.w);
            %obj.w3 =obj.w;
            
        end
        
        function obj = setupTimeAxes(obj)
            obj.t = 0:obj.dt:(obj.n_t-1)*obj.dt;
            obj.n_zp = 2*obj.n_t;
            [obj.T1,obj.T3]=meshgrid(obj.t,obj.t);
        end
        
        function obj = setupSimMatrix(obj)
            obj.simMatrix = zeros(length(obj.w3_in),length(obj.w1_in),length(obj.t2_array)); %array for output
        end
        
        function chi2 = err_fun(obj,p)
            
            obj = obj.calcSpectrum(p);
            chi2 = sum(sum(sum((obj.dataMatrix-obj.simMatrix).^2)));
            
        end
        function chi2 = err_fun_boot(obj,p,ind)
            
            obj = obj.calcSpectrum(p);
            chi2 = sum(sum(sum((obj.dataMatrix(ind)-obj.simMatrix(ind)).^2)));
            
        end
       function chi2 = err_fun_w(obj,p)
            
            obj = obj.calcSpectrum(p);
            chi2 = sum(sum(sum(((obj.dataMatrix-obj.simMatrix).^2).*obj.weightMatrix)));
            
        end
        function chi2 = err_fun_boot_w(obj,p,ind)
            
            obj = obj.calcSpectrum(p);
            chi2 = sum(sum(sum(((obj.dataMatrix(ind)-obj.simMatrix(ind)).^2).*obj.weightMatrix)));
            
        end
        function out = residuals(obj,p)
            if nargin>1
                 obj = obj.calcSpectrum(p);
            end
            out = obj.dataMatrix-obj.simMatrix;
        end
        function out = residuals_w(obj,p)
            if nargin>1
                obj = obj.calcSpectrum(p);
            end
            out = (obj.dataMatrix-obj.simMatrix).*sqrt(obj.weightMatrix);
        end
        function out = residualsSq(obj,p)
            if nargin>1
                obj = obj.calcSpectrum(p);
            end
            out = (obj.dataMatrix-obj.simMatrix).^2;
        end
        function out = residualsSq_w(obj,p)
            if nargin>1
                obj = obj.calcSpectrum(p);
            end
            out = (obj.dataMatrix-obj.simMatrix).^2.*obj.weightMatrix;
        end
        
    end
end
