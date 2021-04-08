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
                obj = obj.validate;
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
        function obj = validate(obj)
            %ensure values of fit parameter are consistent otherwise warn
            if obj.value > obj.ub || obj.value < obj.lb || obj.lb > obj.ub
                warning('SGRLAB:data_analysis:invalid_parameter',...
                    ['fitParamBnd: %s=%f out of bounds (lb = %f, ub = %f) '...
                    'or bounds inconsistent'],...
                    obj.name,obj.value,obj.lb,obj.ub);
            end
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
