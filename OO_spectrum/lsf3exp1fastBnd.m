classdef lsf3exp1fastBnd < lsf3exp1fast & lineshapeFunctionBnd
   properties
       
   end
   
   methods
       function obj = lsf3exp1fastBnd(params,str)
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