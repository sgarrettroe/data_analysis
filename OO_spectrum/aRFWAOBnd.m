classdef aRFWAOBnd < aRFBnd
    
    
    properties
        %
        %  Most properties should be defined in the superclass aRF
        %
        
        %signal processing
        %apodization = 'none'; %### CUT/MOVE
        
        %bandwidth (skip for now) ### CUT/MOVE
        
        % bootstrapping ### CUT/MOVE
        %flag_bootstrap = false;
        %bootstrap_index = [];
        
        %noise (skip for now) ### CUT/MOVE
        
        % orientations ### CUT/MOVE
        %flag_orientational_response = false;
        
        %enumerating diagrams
        rephasing_diagrams;
        nonrephasing_diagrams;
        groundstate_diagrams; %only 01 transitions
        excitedstate_diagrams; %include a 12 transition
        positive_diagrams;
        negative_diagrams;
        
        dyn = additionalDynamics; %structure/class) for additional dynamics
        
        %
        % Additional things that might be fit parameters not in the
        % superclass
        %
        
        
        %transition dipole moments
        mu01sq;
        mu12sqRatio;
        
        %anharmonicity
        anh_cm;
        
    end
    
    methods
        function obj = aRFWAOBnd(options)
            obj@aRFBnd(options);
            %obj = obj.calcSpectrum;
            
            %for subclasses?
            %obj = obj.addNoise;
            %obj = obj.bandwidthFiltering;
            %obj = ...
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

                obj = obj.calcAdditionalDynamics;

                %obj=obj.filterTime;%
                %obj=obj.calcOrientation;
                obj=obj.calcDiagramsFreq(obj.n_zp);
                
                
                %obj=obj.calcAmps;
                %obj=obj.calcKinetics;
                
                %obj=obj.filterFreq;
                
                obj = obj.addDiagrams;
                
                obj=obj.resample(ii);
                %obj=obj.noise;
                
            end
        end
        
        function obj = makeResponseFunctions(obj,~)
            
        end
        function obj = setupResponseFunctions(obj)

            obj.n_diagrams = 8;
            obj.rephasing_diagrams = [1 2 3 7];
            obj.nonrephasing_diagrams = [4 5 6 8];
            obj.groundstate_diagrams = [1 2 4 5 7 8];
            obj.excitedstate_diagrams = [3 6];
            obj.positive_diagrams = [1 2 4 5];
            obj.negative_diagrams = [3 6 7 8];
   
            %have to initialize the first one without an index (don't know
            %why)
            obj.diagrams = feynmanDiagram();
            
            %now calling the last one makes the whole array of empty
            %diagrams
            obj.diagrams(obj.n_diagrams) = feynmanDiagram();

            [obj.diagrams(obj.rephasing_diagrams).isRephasing] = deal(true);
            [obj.diagrams(obj.nonrephasing_diagrams).isRephasing] = deal(false);
            
            
        end
        
        function obj = calcResponseFunctions(obj,p)
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
            
            %update the lineshape function with new parameters
            obj.damping = obj.damping.updateG(p);
            g = obj.damping.g; %shortcut
                        
            rephasingR =  exp(-g(obj.T1)+g(obj.t2)-g(obj.T3)-g(obj.T1+obj.t2)-g(obj.t2+obj.T3)+g(obj.T1+obj.t2+obj.T3));
            nonrephasingR = exp(-g(obj.T1)-g(obj.t2)-g(obj.T3)+g(obj.T1+obj.t2)+g(obj.t2+obj.T3)-g(obj.T1+obj.t2+obj.T3));
            
            [obj.diagrams(obj.rephasing_diagrams).R] = deal(rephasingR);
            [obj.diagrams(obj.nonrephasing_diagrams).R] = deal(nonrephasingR);
            
            ind = obj.negative_diagrams;
            for ii = 1:length(ind)
                obj.diagrams(ind(ii)).R = -obj.diagrams(ind(ii)).R;
            end
        end
        
        function obj = calcPhaseShift(obj,p)

            phi = p.phase_deg/180*pi;
            ind = obj.rephasing_diagrams;
            for ii = 1:length(ind)
                obj.diagrams(ind(ii)).R = ...
                    exp(1i*phi).*obj.diagrams(ind(ii)).R;
            end
            ind = obj.nonrephasing_diagrams;            
            for ii = 1:length(ind)
                obj.diagrams(ind(ii)).R = ...
                    exp(-1i*phi).*obj.diagrams(ind(ii)).R;
            end
            
        end
        
        function obj = calcTDM(obj,p)
            mu_01_2 = p.mu01sq;
            mu_12_2 = p.mu12sqRatio*p.mu01sq;

            phi = p.phase_deg/180*pi;
            for ii = 1:length(obj.groundstate_diagrams)
                obj.diagrams(obj.groundstate_diagrams(ii)).R = ...
                    mu_01_2^2.*obj.diagrams(obj.groundstate_diagrams(ii)).R;
            end
            for ii = 1:length(obj.excitedstate_diagrams)
                obj.diagrams(obj.excitedstate_diagrams(ii)).R = ...
                    mu_01_2*mu_12_2.*exp(1i*phi).*obj.diagrams(obj.excitedstate_diagrams(ii)).R;
            end
        end
        
        function obj = calcAnhShift(obj,p)
            global wavenumbersToInvPs
            anh = p.anh_cm*wavenumbersToInvPs*2*pi;
            term = exp(-1i*anh.*obj.T3);
            ind = obj.excitedstate_diagrams;
            for ii = 1:length(ind)
                obj.diagrams(ind(ii)).R = term.*obj.diagrams(ind(ii)).R;
            end
        end
        
        function obj = calcAdditionalDynamics(obj)
%             persistent T1 T3
% 
%             %only update T1 and T3 if we have to
%             if isempty(T1) || isempty(T3)
%                 T1 = obj.T1;
%                 T3 = obj.T3;
%             end
%             t2 = obj.t2;
            p = obj.paramStruct;
            
            for ii = 1:length(obj.dyn)
                for jj = 1:length(obj.dyn(ii).fun_array)
                    f = obj.dyn(ii).fun_array{jj};
                    ind = obj.dyn(ii).ind_array{jj};
                   for kk = 1:length(ind)
                       obj.diagrams(ind(kk)).R = f(obj.T1,obj.t2,obj.T3,p).*obj.diagrams(ind(kk)).R;
                   end
                end
            end
        end
    end
end

