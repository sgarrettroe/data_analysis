classdef lsfRISDwobbling1NIperp < lineshapeFunction
    
    properties
        params = struct('tr',[],'theta_deg',[]);
        g;
        c2;
        order;
        tpoints;
        L_l;
        R;
    end
    
    methods
        
        function obj = lsfRISDwobbling1NIperp(params,str,aRFoptions)
            if nargin == 0;
                super_args = {};
            else
                super_args{1} = params;
                super_args{2} = str;
            end
            obj@lineshapeFunction(super_args);
        end
        
        function out = makeG(obj)
            
            tr = obj.params(1).tr;
            theta_deg = obj.params(1).theta_deg;
            
            %param struct for R must have these fields tr theta_deg
            p(1).tr = tr;
            p(1).theta_deg = theta_deg;

                        
            %             F_para =@(t, tau) (t-tau).*R.para(p,tau); %this is the FFCF
            F_perp = @(t, tau) (t-tau).*obj.R.perp(p,tau); %this is the FFCF
            
            g_prime = arrayfun(@(t) integral(@(tau) F_perp(t, tau),0,t),obj.tpoints); %do the numerical integration as a function of t
            out = @(t) interp1(obj.tpoints,g_prime,t);
        end
        
        function out = makeC2(obj)
            
            tr = 1/obj.params(1).tr;
            theta_deg = obj.params(1).theta_deg;
            
            p(1).tr = tr;
            p(1).theta_deg = theta_deg;
            
            C = wobblingC;
            Ctot = cell(1,4);
            
            
            for l = 1:4
                %Ctot{l} = 1;
                %for ii = 1:ncones
                Ctot{l}=@(p,t)C{l}(p(1),t);
                %end
            end
            
            out = @(t) 0.4.*(Ctot{2}(p,t2));
            
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

