classdef aRFTwoLevelSystem < aRF
    
    
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
        
        %
        % Additional things that might be fit parameters not in the
        % superclass
        %
                
        
        %transition dipole moments
        mu01sq;
        
        

    end
    
    methods
        function obj = aRFTwoLevelSystem(options)
            obj@aRF(options);
            
            %obj = obj.calcSpectrum;
            
            %for subclasses?
            %obj = obj.addNoise;
            %obj = obj.bandwidthFiltering;
            %obj = ...
        end
               
        function obj = calcSpectrum(obj,p)
            obj = obj.updateFreeFitParams(p);
            obj = obj.setupFreqAxes;
            obj=obj.makeResponseFunctions(obj.paramStruct);
            
            for ii = 1:length(obj.t2_array)
                obj.t2 = obj.t2_array(ii);
                              
                obj=obj.calcDiagramsTime;
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
        
        function obj = makeResponseFunctions(obj,p)
            %update the lineshape function with new parameters
            obj.damping = obj.damping.updateG(p);
            g = obj.damping.g; %shortcut
            
            w_01 = 0; %always 0 in rotating frame%p.w_01_cm*wavenumbersToInvPs*2*pi;
            phi = p.phase_deg/180*pi;
            mu_01_2 = p.mu01sq;
            
            %use
            obj.diagrams = feynmanDiagram(...
                ...%R1(p,l),...
                @(T1,t2,T3)mu_01_2.*exp(-1i*w_01.*(-T1+T3)+1i*phi).*exp(-g(T1)+g(t2)-g(T3)-g(T1+t2)-g(t2+T3)+g(T1+t2+T3)),...
                true);
            
            obj.diagrams(2) = feynmanDiagram(...
                @(T1,t2,T3)mu_01_2.*exp(-1i*w_01.*(T1+T3)-1i*phi).*exp(-g(T1)-g(t2)-g(T3)+g(T1+t2)+g(t2+T3)-g(T1+t2+T3)),...
                false);
            obj.n_diagrams = length(obj.diagrams);
        end
        
        
    end
end

