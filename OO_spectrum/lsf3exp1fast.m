classdef lsf3exp1fast < lineshapeFunction
    
    properties
        params = struct('Delta1_cm',[],'tau1',[],'Delta2_cm',[],'tau2',[],'Delta3_cm',[],'tau3',[],'T2',[]);
        g;
        c2;
    end
    
    methods
        function obj = lsf3exp1fast(params)
            if nargin == 0;
                super_args = {};
            else
                super_args{1} = params;
                super_args{2} = str;
            end
            obj@lineshapeFunction(super_args);
        end
        
        function g = makeG(obj)
            global wavenumbersToInvPs;
        
            Delta1 = obj.params(1).Delta1_cm*wavenumbersToInvPs*2*pi;
            Lambda1 = 1/obj.params(1).tau1;
            Delta2 = obj.params(1).Delta2_cm*wavenumbersToInvPs*2*pi;
            Lambda2 = 1/obj.params(1).tau2;
            Delta3 = obj.params(1).Delta3_cm*wavenumbersToInvPs*2*pi;
            Lambda3 = 1/obj.params(1).tau3;
            T2 = obj.params(1).T2;

            g = @(t) t./T2 + Delta1^2/Lambda1^2.*(exp(-Lambda1.*t)-1+Lambda1.*t)...
                + Delta2^2/Lambda2^2.*(exp(-Lambda2.*t)-1+Lambda2.*t)...
                + Delta3^2/Lambda3^2.*(exp(-Lambda3.*t)-1+Lambda3.*t);
        end
        
        function c2 = makeC2(obj)
            global wavenumbersToInvPs;
        
            Delta1 = obj.params(1).Delta1_cm*wavenumbersToInvPs*2*pi;
            Lambda1 = 1/obj.params(1).tau1;
            Delta2 = obj.params(1).Delta2_cm*wavenumbersToInvPs*2*pi;
            Lambda2 = 1/obj.params(1).tau2;
            Delta3 = obj.params(1).Delta3_cm*wavenumbersToInvPs*2*pi;
            Lambda3 = 1/obj.params(1).tau3;
            T2 = obj.params(1).T2;
            
            c2 = @(t) (t==0)/T2 + Delta1^2.*exp(-Lambda1.*t)...
                + Delta2^2.*exp(-Lambda2.*t)...
                + Delta3^2.*exp(-Lambda3.*t);
        end

    end
end