classdef lsfRISDwobbling1NIperpBnd < lsfRISDwobbling1NIperp & lineshapeFunctionBnd
    properties
        
    end
    
    methods
        function obj = lsfRISDwobbling1NIperpBnd(params,str,aRFoptions)
            if nargin == 0;
                super_args = {};
            else
                super_args{1} = params;
                super_args{2} = str;
                super_args{3} = aRFoptions;
                
            end
            obj@lineshapeFunctionBnd(super_args);
            obj = obj.maketpoints(aRFoptions);
            obj.order = aRFoptions.order;
            obj = obj.makeL_l;
        end
        
        function obj = makeL_l(obj)
            C = wobblingCtest;
            Ctot = cell(1,4);
            
            for l = 1:4
                %Ctot{l} = 1;
                %for ii = 1:ncones
                Ctot{l}=@(p,tau)C{l}(p(1),tau);
                %end
            end
            
            obj.R = wobblingR(Ctot,obj.order);
            obj.L_l = Ctot;
        end        
    end
    
end