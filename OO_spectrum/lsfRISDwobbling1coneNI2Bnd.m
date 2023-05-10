classdef lsfRISDwobbling1coneNI2Bnd < lsfRISDwobbling1coneNI2 & lineshapeFunctionBnd
    properties
        
    end
    
    methods
        function obj = lsfRISDwobbling1coneNI2Bnd(params,str,aRFoptions)
            if nargin == 0
                super_args = {};
            else
                super_args{1} = params;
                super_args{2} = str;
                super_args{3} = aRFoptions;
            end
            obj@lsfRISDwobbling1coneNI2(super_args);
        end
    end
end