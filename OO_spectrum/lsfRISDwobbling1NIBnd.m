classdef lsfRISDwobbling1NIBnd < lsfRISDwobbling1NI & lineshapeFunctionBnd
   properties
       
   end
   
   methods
       function obj = lsfRISDwobbling1NIBnd(params,str,aRFoptions)
            if nargin == 0;
                super_args = {};
            else
                super_args{1} = params;
                super_args{2} = str;
                super_args{3} = aRFoptions;
            end
            obj@lineshapeFunctionBnd(super_args);
            obj = obj.maketpoints(aRFoptions);
        end
 
   end
   
end