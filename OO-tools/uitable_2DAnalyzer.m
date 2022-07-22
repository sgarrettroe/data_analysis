classdef uitable_2DAnalyzer < handle
    
    properties
        UITable

    end

    methods
        function obj = uitable_2DAnalyzer(varargin)

            obj.UITable = uitable(varargin{:});
            
        end

        function rowSelected_Callback(obj)
            
            

        end

    end

end