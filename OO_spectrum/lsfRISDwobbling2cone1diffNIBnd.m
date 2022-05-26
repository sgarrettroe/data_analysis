classdef lsfRISDwobbling2cone1diffNIBnd < lsfRISDwobbling2cone1diffNI & lineshapeFunctionBnd
    properties
        
    end
    
    methods
        function obj = lsfRISDwobbling2cone1diffNIBnd(params,str,aRFoptions)
            if nargin == 0
                super_args = {};
            else
                super_args{1} = params;
                super_args{2} = str;
                super_args{3} = aRFoptions;
            end
            obj@lsfRISDwobbling2cone1diffNI(super_args);
        end
    end
end