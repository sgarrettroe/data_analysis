classdef aDArrayBnd < additionalDynamics & fitParamBnd
    properties
        params; %cell array of parameter structures
        n;
    end
    methods
        function obj = aDArrayBnd(label,fun_array,ind_array,p_array,str)
            if nargin>0
                %could add error checking on inputs here
                obj.label = label;
                obj.n = length(fun_array);
                obj.fun_array = fun_array;
                obj.ind_array = ind_array;
                obj = obj.makeParamStruct(p_array);
                obj = obj.parseStr(str);
            end
        end
        function obj = makeParamStruct(obj,p_array)
            %convert incoming cell array of fitParams arrays to a cell array of structs
            temp = cell2struct({p_array.value},{p_array.name},2);
            temp(2) = cell2struct({p_array.lb},{p_array.name},2);
            temp(3) = cell2struct({p_array.ub},{p_array.name},2);
            obj.params = temp;
        end
        function obj = setParams(obj,pStruct)
            %if the input is a struct
            fnames = fieldnames(pStruct);
            % get names of
            [longnames,shortnames] = regexp(fnames,['(\w+)_' num2str(ii) '$'],'match','tokens');
            ind = cellfun(@(x)~isempty(x),longnames);
            longnames = longnames(ind);
            shortnames = shortnames(ind);
            for jj = 1:length(longnames)
                obj.params(1).(shortnames{jj}{1}{1}) = pStruct.(longnames{jj}{1});
            end
        end
        function obj = setParams2(obj,p)
            %if the input is a fitParam array
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
            for jj = 1:n_fields
                out(jj) = obj.params(1).(fnames{jj});
            end
        end
        function out = paramLb(obj)
            fnames = fieldnames(obj.params);
            n_fields = length(fnames);
            out = zeros(1,length(n_fields));
            for jj = 1:n_fields
                out(jj) = obj.params(2).(fnames{jj});
            end
        end
        function out = paramUb(obj)
            fnames = fieldnames(obj.params);
            n_fields = length(fnames);
            out = zeros(1,length(n_fields));
            for jj = 1:n_fields
                out(jj) = obj.params(3).(fnames{jj});
            end
        end
        
    end
    methods (Access = protected)
        function out = get_name_fxn(obj)
            out = obj.paramNames;
        end
        function out = get_value_fxn(obj)
                out = obj.paramValues;
        end
        function out = get_lb_fxn(obj)
                out = obj.paramLb;
        end
        function out = get_ub_fxn(obj)
            out = obj.paramUb;
        end
        
    end
end