classdef lineshapeFunctionBnd < lineshapeFunction & fitParamBnd
    properties (Hidden)
        lbStruct;
        ubStruct;
    end
    
    methods
        function obj = lineshapeFunctionBnd(args)
            obj@lineshapeFunction({});
            if ~isempty(args)
                params = args{1};
                str = args{2};
                
                obj = obj.setParams2(params);
                obj = obj.parseStr(str);
                
                if ~any( structfun(@isempty, obj.params(1)) )
                    obj.g = makeG(obj);
                    obj.c2 = makeC2(obj);
                end
            end
        end
        function obj = setParams(obj,params)
            fnames = fieldnames(obj.params);
            for ii = 1:length(fnames)
                %set values
                obj.params(1).(fnames{ii}) = params.(fnames{ii});
%                 %set lb as second element of array
%                 obj.params(2).(fnames{ii}) = params.(fnames{ii});                
%                 %set ub as third element of array
%                 obj.params(3).(fnames{ii}) = params.(fnames{ii});

            end
        end
        function obj = setParams2(obj,p)
            fnames = fieldnames(obj.params);
            pnames = {p.name};
            for ii = 1:length(fnames)
                %which elements match by name
                ind = strcmpi(fnames{ii},pnames);
                %set values
                obj.params(1).(fnames{ii}) = p(ind).value;
                %set lb as second element of array
                obj.params(2).(fnames{ii}) = p(ind).lb;                
                %set ub as third element of array
                obj.params(3).(fnames{ii}) = p(ind).ub;
            end
       end
        
    end
        methods (Access = protected)
        function out = get_value_fxn(obj)
            fnames = fieldnames(obj.params);
            n_fields = length(fnames);
            out = zeros(1,length(n_fields));
            for ii = 1:n_fields
                out(ii) = obj.params(1).(fnames{ii});
            end
        end

        function out = get_lb_fxn(obj)
            fnames = fieldnames(obj.params);
            n_fields = length(fnames);
            out = zeros(1,length(n_fields));
            for ii = 1:n_fields
                out(ii) = obj.params(2).(fnames{ii});
            end
        end
        function out = get_ub_fxn(obj)
            fnames = fieldnames(obj.params);
            n_fields = length(fnames);
            out = zeros(1,length(n_fields));
            for ii = 1:n_fields
                out(ii) = obj.params(3).(fnames{ii});
            end
        end
    end

end