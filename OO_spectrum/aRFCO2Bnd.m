classdef aRFCO2Bnd < aRFWAOBnd
    properties
        %
        % Additional values (not fit parameters)
        %
        bend_1q_t0_diagrams;
        bend_1q_t2_diagrams;
        
        %
        % Additional fit params
        %
        dE_cm; %population difference b/n ground and bend excited state
        temperature; %temperature in K
        dw_sb_cm; %stretch-bend coupling counstant in cm-1
        %k_u;
        %k_hgs1;
        %k_hgs2;
    end
    
    methods
        function obj = aRFCO2Bnd(options)
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
                obj = obj.calcBendShift(obj.paramStruct);
                obj = obj.calcThermalPopulations(obj.paramStruct);
                
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
            
            obj.n_diagrams = 32;
            obj.rephasing_diagrams = [1 2 3 7 8 9 13 14 15 19 20 21 25 27 29 31];
            obj.nonrephasing_diagrams = [4 5 6 10 11 12 16 17 18 22 23 24 26 28 30 32];
            obj.groundstate_diagrams = [1 2 4 5 7 8 10 11 13 14 16 17 19 20 22 23 25 26 27 28 29 30 31 32];%all se and gsb peaks
            obj.excitedstate_diagrams = [3 6 9 12 15 18 21 24];
            obj.positive_diagrams = [1 2 4 5 7 8 10 11 13 14 16 17 19 20 22 23]; %blue peaks
            obj.negative_diagrams = [3 6 9 12 15 18 21 24 25:32];%red peaks
            obj.bend_1q_t0_diagrams = [7:12 19:24 29:32];
            obj.bend_1q_t2_diagrams = [7:12 13:18 25 26 29 30];%w_3 = w_bend
            
            %have to initialize the first one without an index (don't know
            %why)
            obj.diagrams = feynmanDiagram();
            
            %now calling the last one makes the whole array of empty
            %diagrams
            obj.diagrams(obj.n_diagrams) = feynmanDiagram();
            
            [obj.diagrams(obj.rephasing_diagrams).isRephasing] = deal(true);
            [obj.diagrams(obj.nonrephasing_diagrams).isRephasing] = deal(false);
            
            
        end
        function obj = calcBendShift(obj,p)
            global wavenumbersToInvPs
            
            dw_sb = p.dw_sb_cm*wavenumbersToInvPs*2*pi;
            termR = exp( 1i*dw_sb.*obj.T1);
            termNR = exp(-1i*dw_sb.*obj.T1);
            ind = obj.bend_1q_t0_diagrams;
            for ii = 1:length(ind)
                if obj.diagrams(ind(ii)).isRephasing
                    obj.diagrams(ind(ii)).R = ...
                        termR.*obj.diagrams(ind(ii)).R;
                else
                    obj.diagrams(ind(ii)).R = ...
                        termNR.*obj.diagrams(ind(ii)).R;
                end
            end
            
            term = exp(-1i*dw_sb.*obj.T3);
            ind = obj.bend_1q_t2_diagrams;
            for ii = 1:length(ind)
                obj.diagrams(ind(ii)).R = ...
                    term.*obj.diagrams(ind(ii)).R;
            end
            
        end
        function obj = calcThermalPopulations(obj,p)
            k_B_cm = 0.69503476; %wavenumbers per K
            p1 = 2*exp(-p.dE_cm/(k_B_cm*p.temperature));
            ind = obj.bend_1q_t0_diagrams;
            for ii = 1:length(ind)
                obj.diagrams(ind(ii)).R = ...
                    p1.*obj.diagrams(ind(ii)).R;
            end
            
        end
        
    end
    
end