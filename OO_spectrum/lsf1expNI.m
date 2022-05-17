classdef lsf1expNI < lineshapeFunction
    
    properties
        params = struct('Delta1_cm',[],'tau1',[]);
        g;
        c2;
    end
    
    methods
        function obj = lsf1expNI(params,str)
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
             F = @(tau,t) (t-tau).*(Delta1^2.*exp(-abs(tau).*Lambda1)); %this is the FFCF
%             F =@(t) Delta1^2/Lambda1^2.*(exp(-Lambda1.*t)-1+Lambda1.*t); %this is the real g(t)
             out =  @(t) arrayfun(@(t) integral(@(tau) F(tau,t),0,t),t); %do the numerical integration as a function of t
        end
        
        function out = makeC2(obj)
            global wavenumbersToInvPs;
        
            Delta1 = obj.params(1).Delta1_cm*wavenumbersToInvPs*2*pi;
            Lambda1 = 1/obj.params(1).tau1;
            
            out = @(t) Delta1^2.*exp(-Lambda1.*t);
        end

    end
end