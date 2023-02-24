classdef lsfRISDwobblingHL1inertial1cone1diff1SSDNI < lineshapeFunction
    
    properties
        params = struct('Delta_cm',[],'theta0_deg',[],'tr1',[],'theta1_deg',[],'tr2',[],'T2',[],'ampSSD',[],'tauSSD',[]);
        g;
        c2;
        order;
        tpoints;
        L_l;
        R;
        pol; %'para' or 'perp'
    end
    
    methods
        
        function obj = lsfRISDwobblingHL1inertial1cone1diff1SSDNI(params,str,aRFoptions) %constructor function
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
                obj = obj.makeL_l;   
            end
        end
        
        function out = makeG(obj)
            global wavenumbersToInvPs;

            % Delta is the total linewidth, not used directly in the RISD
            % response functions, it is in F (scales the whole thing)
            Delta = obj.params(1).Delta_cm*wavenumbersToInvPs*2*pi;
            T2 = obj.params(1).T2;
            a = obj.params(1).ampSSD;
            tauSSD = obj.params(1).tauSSD;            
            %param struct for R must have these fields tr theta_deg
            p = obj.copyParamValuesToParamStruct; 
            
            if strcmpi(obj.pol,'para')
                F =@(t, tau) (t-tau).*a.*exp(-tau./tauSSD).*Delta.^2.*obj.R.para(tau,p); %this is the FFCF time (t-tau) to turn a double integral into a single one
                %F_perp =@(t) (3/25).*(7.*exp(-2.*D_m.*t) - 2.*exp(-12.*D_m.*t)) ./ (1 - 0.4.*exp(-6.*D_m.*t)); %this is the FFCF
            elseif strcmpi(obj.pol,'perp')
                F =@(t, tau) (t-tau).*a.*exp(-tau./tauSSD).*Delta.^2.*obj.R.perp(tau,p); %this is the FFCF time (t-tau) to turn a double integral into a single one
            else
                error('unknown polarization pol = %s, should be either ''para'' or ''perp''\n',obj.pol);
            end
            g_prime = arrayfun(@(t) integral(@(tau) F(t, tau),0,t),obj.tpoints); %do the numerical integration as a function of t
            out = @(t) interp1(obj.tpoints,g_prime,t) + t/T2;
        end
        
        function out = makeC2(obj)
            warning('makeC2 is totally bogus right now and not to be used.')
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
            [C,S] = wobblingCv2;
            Ctot = cell(1,4);            
            for l = 1:4
                %Ctot{l} = 1;
                %for ii = 1:ncones
                Ctot{l}=@(t,p)((S{l}(p.theta0_deg)).^2) ... %inertial cone is just S^2
                    .*C{l}(t,p.tr1,p.theta1_deg)...%reg cone
                    .*exp(-(l*(l+1)./(6.*p.tr2)).*t);%diffusive cone
                %end
            end
            
            obj.R = wobblingRv2(Ctot,obj.order);
            obj.L_l = Ctot;
        end

    end
end