classdef lsf1exp1fastNIBnd < lsf1exp1fastNI & lineshapeFunctionBnd
   properties
       
   end
   
   methods
       function obj = lsf1exp1fastNIBnd(params,str)
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