classdef lsf1exp1fastNI < lineshapeFunction
    
    properties
        params = struct('Delta1_cm',[],'tau1',[],'T2',[]);
        g;
        c2;
        tpoints;
        order;
        pol;
    end
    
    methods
        function obj = lsf1exp1fastNI(params,str,aRFoptions)
            if nargin == 0
                super_args = {};
            elseif nargin == 1 
                super_args = params; %if we were passed cell array
                params = super_args{1};
                str = super_args{2};
                aRFoptions = super_args{3};
            elseif nargin == 2
                super_args{1} = params;
                super_args{2} = str;
                aRFoptions = struct([]);
            elseif nargin == 3
                super_args{1} = params;
                super_args{2} = str;
                super_args{3} = aRFoptions;
            else
                error('confusing number of input args in lsfRISDwobbling1NI: %i\n',nargin)
            end
            obj@lineshapeFunction(super_args);
            if nargin~=0
                %if we have some input arguments
                obj.pol = aRFoptions.pol;
                obj.order = aRFoptions.order;
                obj = obj.maketpoints(aRFoptions);
%                 obj = obj.makeL_l;   
            end           
        end
        
        function g = makeG(obj)
            global wavenumbersToInvPs;
        
            Delta1 = obj.params(1).Delta1_cm*wavenumbersToInvPs*2*pi;
            Lambda1 = 1/obj.params(1).tau1;
            T2 = obj.params(1).T2;
%             g = @(t) t./T2 + Delta1^2/Lambda1^2.*(exp(-Lambda1.*t)-1+Lambda1.*t);
            F = @(tau,t) t./T2 + (t-tau).*(Delta1^2).*exp(-abs(tau).*Lambda1);
            g_prime = arrayfun(@(t) integral(@(tau) F(t, tau),0,t),obj.tpoints); %do the numerical integration as a function of t
            g = @(t) interp1(obj.tpoints,g_prime,t);
        end
        
        function c2 = makeC2(obj)
            global wavenumbersToInvPs;
        
            Delta1 = obj.params(1).Delta1_cm*wavenumbersToInvPs*2*pi;
            Lambda1 = 1/obj.params(1).tau1;
            T2 = obj.params(1).T2;
            
            c2 = @(t) (t==0)/T2 + Delta1^2.*exp(-Lambda1.*t);
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