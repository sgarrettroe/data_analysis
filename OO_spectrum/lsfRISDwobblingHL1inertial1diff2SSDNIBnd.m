classdef lsfRISDwobblingHL1inertial1diff2SSDNIBnd < lsfRISDwobblingHL1inertial1diff2SSDNI & lineshapeFunctionBnd
    properties
        
    end
    
    methods
        function obj = lsfRISDwobblingHL1inertial1diff2SSDNIBnd(params,str,aRFoptions)
            if nargin == 0
                super_args = {};
            else
                super_args{1} = params;
                super_args{2} = str;
                super_args{3} = aRFoptions;
            end
            obj@lsfRISDwobblingHL1inertial1diff2SSDNI(super_args);
        end
    end
end