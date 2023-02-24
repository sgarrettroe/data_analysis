classdef FF3expBnd < FittingFunctionBnd 
   
    properties
        params = struct('amp1',[],'t1',[],'amp2',[],'t2',[],'amp3',[],'t3',[]);
        fit;
    end
    
    methods
           
       function obj = FF3expBnd(params,str)
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
            a2     = obj.params(1).amp2; %variables that are input into the class
            t2     = obj.params(1).t2;
            a3     = obj.params(1).amp3;
            t3     = obj.params(1).t3;
            fit = @(t) a1.*exp(-t./t1) + a2.*exp(-t./t2)+ a3.*exp(-t./t3);
        end
   end
end