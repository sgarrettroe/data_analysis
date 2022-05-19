classdef lsfRISDwobbling1NIpara < lineshapeFunction
    
    properties
        params = struct('Delta1_cm',[],'tr',[],'theta_deg',[]);
        g;
        c2;
        order;
        tpoints;
        L_l;
        R;
        pol='para'; %para or perp
    end
    
    methods
        
        function obj = lsfRISDwobbling1NIpara(params,str,~) %constructor function
            if nargin == 0
                super_args = {};
            else
                super_args{1} = params;
                super_args{2} = str;
            end
            obj@lineshapeFunction(super_args);
        end
        
        function out = makeG(obj)
            tr = obj.params.tr;
            theta_deg = obj.params.theta_deg;
            Delta = obj.Delta;
            
            %param struct for R must have these fields tr theta_deg
            p(1).tr = tr;
            p(1).theta_deg = theta_deg;
            
            if strcmpi(obj.pol,'para')
                F =@(t, tau) (t-tau).*Delta^2.*obj.R.para(p,tau); %this is the FFCF time (t-tau) to turn a double integral into a single one
                %F_perp =@(t) (3/25).*(7.*exp(-2.*D_m.*t) - 2.*exp(-12.*D_m.*t)) ./ (1 - 0.4.*exp(-6.*D_m.*t)); %this is the FFCF
            elseif strcmpi(obj.pol,'perp')
                F =@(t, tau) (t-tau).*Delta^2.*obj.R.perp(p,tau); %this is the FFCF time (t-tau) to turn a double integral into a single one
            else
                error('unknown polarization pol = %s, should be either ''para'' or ''perp''\n',obj.pol);
            end
            g_prime = arrayfun(@(t) integral(@(tau) F(t, tau),0,t),obj.tpoints); %do the numerical integration as a function of t
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

