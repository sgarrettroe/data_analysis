classdef lsfRISDwobbling1NIpara < lineshapeFunction
    
    properties
        params = struct('order',[],'tau_o1',[],'theta_deg1',[]);
        g;
        c2;

        tpoints;
    end
    
    methods
        
        function obj = lsfRISDwobbling1NI(params,str,aRFoptions)
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

                    persistent R old_order;     
                    order = obj.params(1).order;
                    tr1 = 1/obj.params(1).tau_o1;
                    theta_deg1 = obj.params(1).theta_deg1;

                            %param struct for R must have these fields tr theta_deg 
                            p(1).tr = tr1;
                            p(1).theta_deg = theta_deg1;

                            %Look to see if we need to (re)initialize R
                            if isempty(R)||isempty(old_order)||old_order~=order
                                disp('new R');

                                C = wobblingCtest;
                                Ctot = cell(1,4);

                                for l = 1:4
                                    %Ctot{l} = 1;
                                    %for ii = 1:ncones
                                        Ctot{l}=@(p,t)C{l}(p(1),tau);
                                    %end
                                end

                                R = wobblingR(Ctot,order);
                                old_order=order;
                            end
         

            F_para =@(t, tau) (t-tau).*R.para(p,tau); %this is the FFCF
%             F_perp =@(t) (3/25).*(7.*exp(-2.*D_m.*t) - 2.*exp(-12.*D_m.*t)) ./ (1 - 0.4.*exp(-6.*D_m.*t)); %this is the FFCF
           
            g_prime = arrayfun(@(t) integral(@(tau) F_para(t, tau),0,t),obj.tpoints); %do the numerical integration as a function of t
            out = @(t) interp1(obj.tpoints,g_prime,t);
        end
        
        function out = makeC2(obj)
                 persistent Ctot old_order

            order = obj.params(1).order;
            tr1 = 1/obj.params(1).tau_o1;
            theta_deg1 = obj.params(1).theta_deg1;
            
                         p(1).tr = tr1;
                         p(1).theta_deg = theta_deg1;
            %Look to see if we need to (re)initialize R
                if isempty(Ctot)||isempty(old_order)||old_order~=order
                    disp('new Ctot');

                    C = wobblingC;
                    Ctot = cell(1,4);


                    for l = 1:4
                        %Ctot{l} = 1;
                        %for ii = 1:ncones
                            Ctot{l}=@(p,t)C{l}(p(1),t);
                        %end
                    end
                old_order=order;
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

