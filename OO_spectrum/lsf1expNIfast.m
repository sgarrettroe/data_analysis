classdef lsf1expNIfast < lineshapeFunction
    
    properties
        params = struct('Delta1_cm',[],'tau1',[]);
        g;
        c2;

        tpoints;
    end
    
    methods
        function obj = lsf1expNIfast(params,str,aRFoptions)
            if nargin == 0;
                super_args = {};
            else
                super_args{1} = params;
                super_args{2} = str;
            end
            obj@lineshapeFunction(super_args);
            %obj.tpoints = maketpoints(aRFoptions);
        end
        
        function out = makeG(obj)
            global wavenumbersToInvPs;
        
            Delta1 = obj.params(1).Delta1_cm*wavenumbersToInvPs*2*pi;
            Lambda1 = 1/obj.params(1).tau1;
            
            F =@(tau,t) (t-tau).*(Delta1^2.*exp(-abs(tau).*Lambda1)); %this is the FFCF
%             F =@(t) Delta1^2/Lambda1^2.*(exp(-Lambda1.*t)-1+Lambda1.*t); %this is the real g(t)
             %g_prime =  @(t) arrayfun(@(t) integral(@(tau) F(tau,t),0,t),t); %do the numerical integration as a function of t
             g_prime = arrayfun(@(t) integral(@(tau) F(tau,t),0,t),obj.tpoints); %do the numerical integration as a function of t
            out = @(t) interp1(obj.tpoints,g_prime,t);
        end
        
        function out = makeC2(obj)
            global wavenumbersToInvPs;
        
            Delta1 = obj.params(1).Delta1_cm*wavenumbersToInvPs*2*pi;
            Lambda1 = 1/obj.params(1).tau1;
            
            out = @(t) Delta1^2.*exp(-Lambda1.*t);
        end

        function obj = maketpoints(obj,aRFoptions)
            t1 = 0:aRFoptions.dt:(aRFoptions.n_t-1)*aRFoptions.dt;
            t3 = t1;
            t2 = aRFoptions.t2_array;
            tmp = [t1,t3];
            tmp2 = [];
            for ii = 1:length(t2)
                tmp2 = [tmp2,t2(ii), t1 + t2(ii), t2(ii) + t3,t1+t2(ii)+t3];
            end
            obj.tpoints = unique([tmp,tmp2]);
            
        end
    end
end