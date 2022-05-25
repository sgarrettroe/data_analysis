classdef (Abstract) lineshapeFunction < fitParam 
    
    properties (Abstract)
        params;
        g;
        c2;
    end
    
    methods (Abstract)
        makeG(obj);
        makeC2(obj);
    end
    
    methods
        function obj = lineshapeFunction(args)
            %            obj@fitParam(args{1});
            %            if length(args)>1
            %                params = args{2};
            
            %obj.name = class(obj);
            
            if ~isempty(args)
                params = args{1};
                str = args{2};
                
                obj = obj.parseStr(str);
                if isa(params,'fitParam')
                    obj = obj.setParams2(params);                    
%                     if ~any(isempty(obj.params))
%                         obj.g = makeG(obj);
%                         obj.c2 = makeC2(obj);
%                     end
                else
                    obj = obj.setParams(params);
%                     if ~any( structfun(@isempty, obj.params) )
%                         obj.g = makeG(obj);
%                         obj.c2 = makeC2(obj);
%                     end
                end
                
            end
        end
        
        function obj = setParams(obj,params)
            fnames = fieldnames(obj.params);
            for ii = 1:length(fnames)
                obj.params.(fnames{ii}) = params.(fnames{ii});
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

        function out = paramNames(obj)
            %f = @(nameCell) horzcat('damping.',nameCell);
            out = fieldnames(obj.params)';
            %             for ii = 1:length(out),
            %                 out{ii} = horzcat('damping.',out{ii});
            %             end
        end
        
        function out = paramValues(obj)
            fnames = fieldnames(obj.params);
            n_fields = length(fnames);
            out = zeros(1,length(n_fields));
            for ii = 1:n_fields
                out(ii) = obj.params.(fnames{ii});
            end
        end
        
        function obj = updateG(obj,p)
            obj = obj.setParams(p);
            obj.g = obj.makeG;
        end
    end
    methods (Access = protected)
        function out = get_name_fxn(obj)
            out = obj.paramNames;
        end
        function out = get_value_fxn(obj)
            out = obj.paramValues;
        end
    end
end