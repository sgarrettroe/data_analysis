classdef FF2VoigtBnd < FittingFunctionBnd
   
    properties
        params = struct('amp1',[],'amp2',[],'w_gauss',[],'w_lorentz',[],'center',[],'anh',[]);
        fit;
    end
    
       methods
           
       function obj = FF2VoigtBnd(params,str)
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
            a2     = obj.params(1).amp2;
            wg     = obj.params(1).w_gauss;
            wl     = obj.params(1).w_lorentz;
            center = obj.params(1).center;
            anh    = obj.params(1).anh;
            
            fit = @(w) - a1.*voigt(w,center,wg,wl) + a2.*voigt(w,center - anh,wg,wl);
        end
   end
end