classdef lsf1expBndNI < lsf1expNI & lineshapeFunctionBnd
   properties
       
   end
   
   methods
       function obj = lsf1expBndNI(params,str)
            if nargin == 0;
                super_args = {};
            else
                super_args{1} = params;
                super_args{2} = str;
            end
            obj@lineshapeFunctionBnd(super_args);
        end
 
   end
   
end