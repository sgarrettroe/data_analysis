classdef lsfRISDwobbling1inertial1diffNIBnd < lsfRISDwobbling1inertial1diffNI & lineshapeFunctionBnd
    properties
        
    end
    
    methods
        function obj = lsfRISDwobbling1inertial1diffNIBnd(params,str,aRFoptions)
            if nargin == 0
                super_args = {};
            else
                super_args{1} = params;
                super_args{2} = str;
                super_args{3} = aRFoptions;
            end
            obj@lsfRISDwobbling1inertial1diffNI(super_args);
        end
    end
end