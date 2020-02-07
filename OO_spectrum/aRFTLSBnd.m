classdef aRFTLSBnd < aRFTwoLevelSystem & aRFBnd
    properties

    end
    methods
        function obj = aRFTLSBnd(options)
            obj@aRFBnd(options);
            obj@aRFTwoLevelSystem(options);
        end
    end
    
end
