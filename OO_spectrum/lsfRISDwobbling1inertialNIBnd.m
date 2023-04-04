classdef lsfRISDwobbling1inertialNIBnd < lsfRISDwobbling1inertialNI & lineshapeFunctionBnd
    properties
        
    end
    
    methods
        function obj = lsfRISDwobbling1inertialNIBnd(params,str,aRFoptions)
            if nargin == 0
                super_args = {};
            else
                super_args{1} = params;
                super_args{2} = str;
                super_args{3} = aRFoptions;
            end
            obj@lsfRISDwobbling1inertialNI(super_args);
        end
    end
end