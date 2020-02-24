classdef fitParam 
    properties
        name;
        value;
        isFree;
    end
    properties (Hidden)
        nameHidden;
        valueHidden;
    end
    methods
        function obj = fitParam(name,val,stringFreeOrFixed)
            if nargin>0
                obj.name = name;
                obj.value = val;
                obj = obj.parseStr(stringFreeOrFixed);
                
            end
            %             if nargin>1
            %                 obj.value = val;
            %             end
        end
        function out = get.name(obj)
            out = obj.get_name_fxn;
        end
        function out = get.value(obj)
            out = obj.get_value_fxn;
        end
        function obj = set.name(obj,val)
            obj = obj.set_name_fxn(val);
        end
        function obj = set.value(obj,val)
            obj = obj.set_value_fxn(val);
        end
        
        function obj = parseStr(obj,str)
                switch lower(str)
                    case {'free'}
                        obj.isFree = true;
                    case {'fixed'}
                        obj.isFree = false;
                    case {''}
                        obj.isFree = [];
                    otherwise
                        error('parameter must be either free or fixed (or empty for grouped parameters like lineshapeFunctions)');
                end
        end            
    end
    methods (Access = protected) %these can be customized
        function out = get_name_fxn(obj)
            out = obj.nameHidden;
        end
        function out = get_value_fxn(obj)
            out = obj.valueHidden;
        end
        function obj = set_name_fxn(obj,val)
            obj.nameHidden = val;
        end
        function obj = set_value_fxn(obj,val)
            obj.valueHidden = val;
        end
    end
end