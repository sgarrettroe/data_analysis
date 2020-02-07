classdef lsf1exp < lineshapeFunction
    
    properties
        params = struct('Delta1_cm',[],'tau1',[]);
        g;
        c2;
    end
    
    methods
        function obj = lsf1exp(params,str)
            if nargin == 0;
                super_args = {};
            else
                super_args{1} = params;
                super_args{2} = str;
            end
            obj@lineshapeFunction(super_args);
        end
        
        function out = makeG(obj)
            global wavenumbersToInvPs;
        
            Delta1 = obj.params(1).Delta1_cm*wavenumbersToInvPs*2*pi;
            Lambda1 = 1/obj.params(1).tau1;

            out = @(t) Delta1^2/Lambda1^2.*(exp(-Lambda1.*t)-1+Lambda1.*t);
        end
        
        function out = makeC2(obj)
            global wavenumbersToInvPs;
        
            Delta1 = obj.params(1).Delta1_cm*wavenumbersToInvPs*2*pi;
            Lambda1 = 1/obj.params(1).tau1;
            
            out = @(t) Delta1^2.*exp(-Lambda1.*t);
        end

    end
end