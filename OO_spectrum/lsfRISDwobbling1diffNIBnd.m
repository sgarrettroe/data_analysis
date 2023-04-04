classdef lsfRISDwobbling1diffNIBnd < lsfRISDwobbling1diffNI & lineshapeFunctionBnd
    properties
        
    end
    
    methods
        function obj = lsfRISDwobbling1diffNIBnd(params,str,aRFoptions)
            if nargin == 0
                super_args = {};
            else
                super_args{1} = params;
                super_args{2} = str;
                super_args{3} = aRFoptions;
            end
            obj@lsfRISDwobbling1diffNI(super_args);
        end
    end
end