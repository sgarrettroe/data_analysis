classdef lsfRISDwobblingFieldHL1cone1diff1SSDNIBnd < lsfRISDwobblingFieldHL1cone1diff1SSDNI & lineshapeFunctionBnd
    properties
        
    end
    
    methods
        function obj = lsfRISDwobblingFieldHL1cone1diff1SSDNIBnd(params,str,aRFoptions)
            if nargin == 0
                super_args = {};
            else
                super_args{1} = params;
                super_args{2} = str;
                super_args{3} = aRFoptions;
            end
            obj@lsfRISDwobblingFieldHL1cone1diff1SSDNI(super_args);
        end
    end
end