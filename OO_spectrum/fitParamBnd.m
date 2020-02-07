classdef fitParamBnd < fitParam
    
    properties
        lb;
        ub;
    end
    properties (Hidden)
        lbHidden;
        ubHidden;
    end

    methods
        function obj = fitParamBnd(name,value,lb,ub,freeOrFixed)
            if nargin>0
                obj.name = name;
                obj.value = value;
                obj.lb = lb;
                obj.ub = ub;
                obj = obj.parseStr(freeOrFixed);
            end
        end
        function out = get.lb(obj)
            out = obj.get_lb_fxn;
        end
        function out = get.ub(obj)
            out = obj.get_ub_fxn;
        end
        function obj = set.lb(obj,val)
            obj = obj.set_lb_fxn(val);
        end
        function obj = set.ub(obj,val)
            obj = obj.set_ub_fxn(val);
        end
        
    end
    methods (Access = protected) %these can be customized
        function out = get_lb_fxn(obj)
            out = obj.lbHidden;
        end
        function out = get_ub_fxn(obj)
            out = obj.ubHidden;
        end
        function obj = set_lb_fxn(obj,val)
            obj.lbHidden = val;
        end
        function obj = set_ub_fxn(obj,val)
            obj.ubHidden = val;
        end
    end
    
end
