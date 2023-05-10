classdef lsfRISDwobblingHL1cone1SSDNIBnd < lsfRISDwobblingHL1cone1SSDNI & lineshapeFunctionBnd
    properties
        
    end
    
    methods
        function obj = lsfRISDwobblingHL1cone1SSDNIBnd(params,str,aRFoptions)
            if nargin == 0
                super_args = {};
            else
                super_args{1} = params;
                super_args{2} = str;
                super_args{3} = aRFoptions;
            end
            obj@lsfRISDwobblingHL1cone1SSDNI(super_args);
        end
    end
end