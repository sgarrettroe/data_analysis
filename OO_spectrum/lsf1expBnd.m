classdef lsf1expBnd < lsf1exp & lineshapeFunctionBnd
   properties
       
   end
   
   methods
       function obj = lsf1expBnd(params,str)
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