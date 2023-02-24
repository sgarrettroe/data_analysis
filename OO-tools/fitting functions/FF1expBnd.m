classdef FF1expBnd < FittingFunctionBnd 
   
    properties
        params = struct('amp1',[],'t1',[]);
        fit;
    end
    
    methods
           
       function obj = FF1expBnd(params,str)
            if nargin == 0;
                super_args = {};
            else
                super_args{1} = params;
                super_args{2} = str;
            end
            obj@FittingFunctionBnd(super_args);
       end
       
       function fit = makeFit(obj) %makes the appropriate fitting function
        
            a1     = obj.params(1).amp1; %variables that are input into the class
            t1     = obj.params(1).t1;

            fit = @(t) a1.*exp(-t./t1);
        end
   end
end