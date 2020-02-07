classdef lsfArrayBnd < lineshapeFunctionBnd
    properties
        params;
        g;
        c2;
        lsf_array; %cell array of lsfs
        ind_array; %list of what diagrams the lineshape function applies to
        n; %number of elements
    end
    methods
        function obj = lsfArrayBnd(lsf_cell_in,ind_cell_in,string)
            obj@lineshapeFunctionBnd({});
            if nargin>0
                obj.lsf_array = lsf_cell_in;
                obj.ind_array = ind_cell_in;
                obj.n = length(obj.lsf_array);
                obj = obj.parseStr(string);
            end
        end
        function obj = makeG(obj)
            for ii=1:obj.n
                obj.lsf_array{ii} =  obj.lsf_array{ii}.makeG;
            end
        end
        function obj = makeC2(obj)
            for ii=1:obj.n
                obj.lsf_array{ii} =  obj.lsf_array{ii}.makeC2;
            end
        end
        function obj = updateG(obj,pStruct)
            pCell = separateParams(obj,pStruct);
            for ii = 1:obj.n
                obj.lsf_array{ii} =  obj.lsf_array{ii}.updateG(pCell{ii});
            end
        end
        
        function out = separateParams(obj,pStruct)
            fnames = fieldnames(pStruct);
            out = cell(1,obj.n);
            for ii =1:obj.n
                % get names of
                [longnames,shortnames] = regexp(fnames,['(\w+)_' num2str(ii) '$'],'match','tokens');
                ind = cellfun(@(x)~isempty(x),longnames);
                longnames = longnames(ind);
                shortnames = shortnames(ind);
                for jj = 1:length(longnames)
                    out{ii}.(shortnames{jj}{1}{1}) = pStruct.(longnames{jj}{1});
                end
            end
        end
    end
    methods (Access = protected)
        function out = get_name_fxn(obj)
            out = [];
            for ii = 1:obj.n
                str = ['_' num2str(ii)];
                out = [out strcat(obj.lsf_array{ii}.paramNames,str)];
            end
        end
        function out = get_value_fxn(obj)
            out = [];
            for ii = 1:obj.n
                out = [out obj.lsf_array{ii}.paramValues];
            end
        end
        function out = get_lb_fxn(obj)
            out = [];
            for ii =1:obj.n
                out = [out obj.lsf_array{ii}.lb];
            end
        end
        function out = get_ub_fxn(obj)
            out = [];
            for ii =1:obj.n
                out = [out obj.lsf_array{ii}.ub];
            end
        end
        
    end
end


%
%    properties (Abstract)
%         params;
%         g;
%         c2;
%     end
%
%     methods (Abstract)
%         makeG(obj);
%         makeC2(obj);
%     end
%