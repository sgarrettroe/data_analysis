classdef additionalDynamics
    properties
        label;
        fun_array;
        ind_array;
    end
    methods
        function obj = additionalDynamics(label,fun_array,ind_array)
            if nargin>0
                %could add error checking on inputs here
                obj.label = label;
                obj.fun_array = fun_array;
                obj.ind_array = ind_array;
            end
        end
    end
end