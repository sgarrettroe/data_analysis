classdef feynmanDiagram
    
    properties
        fun; %for function handle
        isRephasing;
        R; %calculated response
    end
    
    methods
        function obj = feynmanDiagram(fun,isRephasing)
            if nargin>0
                obj.fun = fun;
                obj.isRephasing = isRephasing;
            end
        end
        
        function obj = calcResponseTime(obj,T1,t2,T3)
            obj.R = obj.fun(T1,t2,T3);
        end
        
        function obj = timeToFreq(obj,n_zp)
            if obj.isRephasing
                obj = timeToFreqRephasing(obj,n_zp);
            else
                obj = timeToFreqNonRephasing(obj,n_zp);
            end
        end
        
        function obj = timeToFreqRephasing(obj,n_zp)
            
            obj.R = sgrsfft2(obj.R,n_zp);
            obj.R = -fftshift(real(fliplr(circshift(obj.R,[0 -1]))));
            
        end
        
        function obj = timeToFreqNonRephasing(obj,n_zp)
            
            obj.R = sgrsfft2(obj.R,n_zp);
            obj.R = -fftshift(real(obj.R));
            
        end
    end
    
end
