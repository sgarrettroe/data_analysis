classdef lsfRISDwobbling1NI < lineshapeFunction
    
    properties
        params = struct('Delta_cm',[],'tr1',[],'theta1_deg',[]);
        g;
        c2;
        order;
        tpoints;
        L_l;
        R;
        pol; %'para' or 'perp'
    end
    
    methods
        
        function obj = lsfRISDwobbling1NI(params,str,aRFoptions) %constructor function
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
            elseif naragin == 3
                super_args{1} = params;
                super_args{2} = str;
                super_args{3} = aRFoptions;
            else
                error('confusing number of input args in lsfRISDwobbling1NIpara: %i\n',nargin)
            end
            obj@lineshapeFunction(super_args);
            if nargin~=0
                %if we have some input arguments
                obj.pol = aRFoptions.pol;
                obj.order = aRFoptions.order;
                obj = obj.maketpoints(aRFoptions);
                obj = obj.makeL_l;   
            end
        end
        
        function out = makeG(obj)
            global wavenumbersToInvPs;

            tr1 = obj.params(1).tr1;
            theta1_deg = obj.params(1).theta1_deg;
            Delta = obj.params(1).Delta_cm*wavenumbersToInvPs*2*pi;
            
            %param struct for R must have these fields tr theta_deg
            p(1).tr1 = tr1;
            p(1).theta1_deg = theta1_deg;
            
            if strcmpi(obj.pol,'para')
                F =@(t, tau) (t-tau).*Delta^2.*obj.R.para(tau,p); %this is the FFCF time (t-tau) to turn a double integral into a single one
                %F_perp =@(t) (3/25).*(7.*exp(-2.*D_m.*t) - 2.*exp(-12.*D_m.*t)) ./ (1 - 0.4.*exp(-6.*D_m.*t)); %this is the FFCF
            elseif strcmpi(obj.pol,'perp')
                F =@(t, tau) (t-tau).*Delta^2.*obj.R.perp(tau,p); %this is the FFCF time (t-tau) to turn a double integral into a single one
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
            
            C = wobblingCv2;
            Ctot = cell(1,4);
            
            for l = 1:4
                %Ctot{l} = 1;
                %for ii = 1:ncones
                Ctot{l}=@(p,t)C{l}(t,p(1));
                %end
            end
            
            out = @(t) 0.4.*(Ctot{2}(t,p));
            
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
        function obj = makeL_l(obj)
            C = wobblingCv2;
            Ctot = cell(1,4);            
            for l = 1:4
                %Ctot{l} = 1;
                %for ii = 1:ncones
                Ctot{l}=@(tau,p)C{l}(tau,p.tr1,p.theta1_deg);
                %end
            end
            
            obj.R = wobblingRv2(Ctot,obj.order);
            obj.L_l = Ctot;
        end

    end
end

