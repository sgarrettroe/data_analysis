classdef aRFEANBnd < aRFWAOBnd
    properties
        %
        % Additional values (not fit parameters)
        %
        low_freq_diagrams;
        
        %
        % Additional fit params
        %
        amp2; %amplitude of low freq peak
        dw_cm; %shift of low freq peak in cm-1
        tau_o_1;
        tau_o_2;
        T1_1;
        T1_2;
    end
    
    methods
        function obj = aRFEANBnd(options)
            obj@aRFWAOBnd(options);
            
        end
        function obj = calcSpectrum(obj,p)
            obj = obj.updateFreeFitParams(p);
            obj = obj.setupFreqAxes;
            %obj=obj.makeResponseFunctions(obj.paramStruct);
            obj=obj.setupResponseFunctions;
            
            for ii = 1:length(obj.t2_array)
                obj.t2 = obj.t2_array(ii);
                
                %obj=obj.calcDiagramsTime;
                obj = obj.calcResponseFunctions(obj.paramStruct);%could remove input parameter...
                
                obj = obj.calcPhaseShift(obj.paramStruct);
                obj = obj.calcAnhShift(obj.paramStruct);
                obj = obj.calcTDM(obj.paramStruct);
                obj = obj.calcFreqShift(obj.paramStruct);
                obj = obj.calcAmp(obj.paramStruct);
                
                obj = obj.calcAdditionalDynamics;
                
                %obj=obj.filterTime;%
                obj=obj.calcDiagramsFreq(obj.n_zp);
                
                %obj=obj.filterFreq;
                
                obj = obj.addDiagrams;
                
                obj=obj.resample(ii);
                %obj=obj.noise;
                
            end
        end
        
        function obj = makeResponseFunctions(obj,~)
            
        end
        function obj = setupResponseFunctions(obj)
            
            obj.n_diagrams = 12;
            obj.rephasing_diagrams = [1 2 3 7 8 9];
            obj.nonrephasing_diagrams = [4 5 6 10 11 12];
            obj.groundstate_diagrams = [1 2 4 5 7 8 10 11];
            obj.excitedstate_diagrams = [3 6 9 12];
            obj.positive_diagrams = [1 2 4 5 7 8 10 11];
            obj.negative_diagrams = [3 6 9 12];
            obj.low_freq_diagrams = [7:12];
            
            %have to initialize the first one without an index (don't know
            %why)
            obj.diagrams = feynmanDiagram();
            
            %now calling the last one makes the whole array of empty
            %diagrams
            obj.diagrams(obj.n_diagrams) = feynmanDiagram();
            
            [obj.diagrams(obj.rephasing_diagrams).isRephasing] = deal(true);
            [obj.diagrams(obj.nonrephasing_diagrams).isRephasing] = deal(false);
            
            
        end
        function obj = calcFreqShift(obj,p)
            global wavenumbersToInvPs
            
            dw = p.dw_cm*wavenumbersToInvPs*2*pi;
            termNR = exp(-1i*dw.*obj.T1).*exp(-1i*dw.*obj.T3);
            termR  = exp( 1i*dw.*obj.T1).*exp(-1i*dw.*obj.T3);
            ind = obj.low_freq_diagrams;
            for ii = 1:length(ind)
                if obj.diagrams(ind(ii)).isRephasing
                    obj.diagrams(ind(ii)).R = ...
                        termR.*obj.diagrams(ind(ii)).R;
                else
                    obj.diagrams(ind(ii)).R = ...
                        termNR.*obj.diagrams(ind(ii)).R;
                end
            end
            
        end
        function obj = calcAmp(obj,p)
            amp = p.amp2;
            ind = obj.low_freq_diagrams;
            for ii = 1:length(ind)
                obj.diagrams(ind(ii)).R = ...
                    amp.*obj.diagrams(ind(ii)).R;
            end
        end
        function obj = calcResponseFunctions(obj,p)
            % modified to take array of inputs
            
            %             persistent T1 T3
            %
            %             %only update T1 and T3 if we have to
            %             if isempty(T1) || isempty(T3)
            %                 T1 = obj.T1;
            %                 T3 = obj.T3;
            %             end
            %
            %             %always update t2
            %             t2 = obj.t2;
            
            if isa(obj.damping,'lsfArrayBnd')
                for ii = 1:obj.damping.n
                    
                    %update the lineshape function with new parameters
                    obj.damping = obj.damping.updateG(p);
                    
                    g = obj.damping.lsf_array{ii}.g; %shortcut
                    
                    rephasingR =  exp(-g(obj.T1)+g(obj.t2)-g(obj.T3)-g(obj.T1+obj.t2)-g(obj.t2+obj.T3)+g(obj.T1+obj.t2+obj.T3));
                    nonrephasingR = exp(-g(obj.T1)-g(obj.t2)-g(obj.T3)+g(obj.T1+obj.t2)+g(obj.t2+obj.T3)-g(obj.T1+obj.t2+obj.T3));
                    
                    ind = obj.damping.ind_array{ii};
                    for jj = 1:length(ind)
                        if obj.diagrams(ind(jj)).isRephasing
                            %                    [obj.diagrams(obj.rephasing_diagrams).R] = deal(rephasingR);
                            obj.diagrams(ind(jj)).R = rephasingR;
                        else
                            %                   [obj.diagrams(obj.nonrephasing_diagrams).R] = deal(nonrephasingR);
                            obj.diagrams(ind(jj)).R = nonrephasingR;
                        end
                        
                    end
                    
                    
                end
                ind = obj.negative_diagrams;
                for ii = 1:length(ind)
                    obj.diagrams(ind(ii)).R = -obj.diagrams(ind(ii)).R;
                end
            else
                % if we have the usual input, call the usual function
                obj = calcResponseFunctions@aRFWAOBnd(obj,p);
            end
        end
    end
end