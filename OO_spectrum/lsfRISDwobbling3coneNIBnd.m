classdef lsfRISDwobbling3coneNIBnd < lsfRISDwobbling3coneNI & lineshapeFunctionBnd
    properties
        
    end
    
    methods
        function obj = lsfRISDwobbling3coneNIBnd(params,str,aRFoptions)
            if nargin == 0
                super_args = {};
            else
                super_args{1} = params;
                super_args{2} = str;
                super_args{3} = aRFoptions;
            end
            obj@lsfRISDwobbling3coneNI(super_args);
        end
    end
end