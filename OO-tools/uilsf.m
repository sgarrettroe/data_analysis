classdef uilsf < handle
    
    properties (Access = public)
        Parent
        DropDown matlab.ui.control.DropDown
        DropDownLabel matlab.ui.control.Label %deprecated
        OuterGrid matlab.ui.container.GridLayout
        InnerGrid matlab.ui.container.GridLayout
        LoadProfileButton
        SaveProfileButton
        DoneButton
        %lsfProfile Profile %deprecated
        lsf
        s = struct('vals',[],'lb',[],'ub',[],'isFree',[]);
    end
    
    properties (SetAccess = immutable)
        PROPERTIES_TO_SAVE = {'lsf'}
    end
    
    methods (Access = public)
        function obj = uilsf(Parent)
            %function obj = uilsf(Parent,Grid)
            %function obj = uilsf(Grid)
            obj.Parent = Parent;
            %obj.OuterGrid = Grid;
            obj = createComponents(obj);
        end
        function updateLsfDropDown(obj,dropDown)
            searchTarget = 'lsf1exp'; %a file that should exist to help us find the folder
            if exist(searchTarget,'file')
                path = dir(which(searchTarget)).folder;
                lsfNames = {dir([path,filesep,'lsf*Bnd.m']).name};
            end
            dropDown.Items = [{'Select One'},lsfNames(:)'];
        end
        
        %         function updateLsfTable(obj)
        %             disp('updating lsf table')
        %             if exist(obj.lsfName,'file')
        %                 obj.lsf = feval(obj.lsfName);
        %                 %lsfTable.RowName = obj.lsf.name;
        %                 col1 = obj.lsf.name';ed
        %                 %                 col2 = double.empty(length(obj.lsf.name),0); %sp
        %                 %                 col3 = double.empty(length(obj.lsf.name),0); %lb
        %                 %                 col4 = double.empty(length(obj.lsf.name),0); %ub
        %                 col2 = nan(length(obj.lsf.name),1); %sp
        %                 col3 = nan(length(obj.lsf.name),1); %lb
        %                 col4 = nan(length(obj.lsf.name),1); %ub
        %                 t = table(col1,col2,col3,col4);
        %                 lsfTable.ColumnEditable=[false,true,true,true];
        %                 lsfTable.Data = t;
        %             else
        %                 warning(['lineshape function ',obj.lsfName,' not found.']);
        %             end
        %         end
        
        function updateLsfInnerGrid(obj)
            disp('updateLsfInnerGrid');
            
            % get all the parameter names in cell array
            names = obj.lsf.name;
            try
                % if it isn't empty this should work
                vals = obj.lsf.value;
                lb = obj.lsf.lb;
                ub = obj.lsf.ub;
                if strcmpi(obj.lsf.isFree,'free')
                    isFree = true;
                else
                    isFree = false;
                end
            catch
                % if it was reset then set all to 0
                %warning(E);
                vals = zeros(1,length(names));
                lb = zeros(1,length(names));
                ub = zeros(1,length(names));
                isFree = false;
            end
            obj.s.vals = vals;
            obj.s.lb = lb;
            obj.s.ub = ub;
            obj.s.isFree = isFree;
            
            %clear previous entries
            g = obj.InnerGrid;
            if ~isempty(g),delete(g.Children);end
            
            % how many params <=> nrows
            nnames = length(names);
            nrows = nnames+2; %+ headers and isFree rows
            
            % set grid rows and columns
            g.RowHeight = repmat({22},1,nrows);
            
            %
            heads = {'','sp','lb','ub'};
            ncols = length(heads); %label, sp, lb, ub
            g.ColumnWidth = [{100},repmat({'1x'},1,ncols-1)];
            
            % make the headings
            columnHeadings = cell(1,ncols);
            for jj = 1:ncols
                columnHeadings{jj}=uilabel(g);
                columnHeadings{jj}.Layout.Row = 1;
                columnHeadings{jj}.Layout.Column = jj;
                columnHeadings{jj}.Text = heads{jj};
                columnHeadings{jj}.HorizontalAlignment = 'center';
            end
            
            labels = cell(nrows,1);
            editFields = cell(nrows,ncols);
            for ii = 1:nnames
                %labels
                labels{ii} = uilabel(g);
                labels{ii}.Layout.Row = ii+1;
                labels{ii}.Layout.Column = 1;
                labels{ii}.Text = names{ii};
                labels{ii}.HorizontalAlignment = 'right';
                
                %sp
                editFields{ii,1} = uieditfield(g,'numeric');
                editFields{ii,1}.Layout.Row = ii+1;
                editFields{ii,1}.Layout.Column = 2;
                editFields{ii,1}.Value = vals(ii);
                editFields{ii,1}.ValueChangedFcn = @(src,event) EditFieldValueChanged(obj,src,event,ii,1);
                
                %lb
                editFields{ii,2} = uieditfield(g,'numeric');
                editFields{ii,2}.Layout.Row = ii+1;
                editFields{ii,2}.Layout.Column = 3;
                editFields{ii,2}.Value = lb(ii);
                editFields{ii,2}.ValueChangedFcn = @(src,event) EditFieldValueChanged(obj,src,event,ii,2);
                
                %ub
                editFields{ii,3} = uieditfield(g,'numeric');
                editFields{ii,3}.Layout.Row = ii+1;
                editFields{ii,3}.Layout.Column = 4;
                editFields{ii,3}.Value = ub(ii);
                editFields{ii,3}.ValueChangedFcn = @(src,event) EditFieldValueChanged(obj,src,event,ii,3);
                
            end
            
            chkbox = uicheckbox(g);
            chkbox.Text = 'is free';
            chkbox.Value = isFree;
            chkbox.Layout.Row = nrows; %bottom
            chkbox.Layout.Column = 2; %left?
            chkbox.ValueChangedFcn = @(src,event) CheckBoxValueChanged(obj,src,event);
            
            %lsfTable.RowName = obj.lsf.name;
            %                 col1 = obj.lsf.name';
            %                 %                 col2 = double.empty(length(obj.lsf.name),0); %sp
            %                 %                 col3 = double.empty(length(obj.lsf.name),0); %lb
            %                 %                 col4 = double.empty(length(obj.lsf.name),0); %ub
            %                 col2 = nan(length(obj.lsf.name),1); %sp
            %                 col3 = nan(length(obj.lsf.name),1); %lb
            %                 col4 = nan(length(obj.lsf.name),1); %ub
            %                 t = table(col1,col2,col3,col4);
            %                 lsfTable.ColumnEditable=[false,true,true,true];
            %                 lsfTable.Data = t;
        end
        
        
        % Drop down opening function: DropDown
        function DropDownOpening(obj, src)
            searchTarget = 'lsf1exp'; %a file that should exist to help us find the folder
            if exist(searchTarget,'file')
                path = dir(which(searchTarget)).folder;
                lsfNames = {dir([path,filesep,'lsf*Bnd.m']).name};
            end
            src.Items = [{'Select One'},lsfNames(:)'];
        end
        
        % Value changed function: DropDown
        function DropDownValueChanged(obj, src)
            %clear old lsf
            obj.lsf = [];
            
            % build empty new one
            [~,lsfName,~] = fileparts([filesep,src.Value]);
            
            if exist(lsfName,'file')
                obj.lsf = feval(lsfName);
            else
                warning(['SGRLAB:lineshape function ',lsfName,' not found.']);
            end
            
            % update grid accordingly
            updateLsfInnerGrid(obj);
        end
        
        function EditFieldValueChanged(obj,src,event,ii,jj)
            switch jj
                case 1
                    obj.s.vals(ii) = event.Value;
                case 2
                    obj.s.lb(ii) = event.Value;
                case 3
                    obj.s.ub(ii) = event.Value;
                otherwise
                    warning('SGRLAB: value of jj in EditFieldValueChanged out of bounds.')
            end
        end
        
        function CheckBoxValueChanged(obj,src,event)
            obj.s.isFree = event.Value;
        end
        
        function updateLsf(obj)
            % take the values out of s.vals, s.lb, s.ub (which come from
            % the `table', and save them in the lsf obj
            names = obj.lsf.name;
            nparams = length(names);
            for ii = nparams:-1:1 %going backwards so it is preallocated
                lsfParams(ii) = fitParamBnd(names{ii},obj.s.vals(ii),obj.s.lb(ii),obj.s.ub(ii),'');
            end
            obj.lsf = obj.lsf.setParams2(lsfParams);
            if obj.s.isFree
                obj.lsf.isFree = 'free';
            else
                obj.lsf.isFree = 'fixed';
            end
        end
        
        function loadProfile(obj,src,event)
            %load a profile of the class (or a class) by name
            
            [filename,pathname] = uigetfile('*.mat','Select a file to load');
            if filename==0
                return
            end
            file_and_path = [pathname, filesep,filename];
            if exist(file_and_path,'file')
                tmps = load(file_and_path);
            else
                warning(['SGRLAB:IO','Could not find file ', [pathname, filesep,filename],' to load.']);
            end
            
            thesefieldnames = fieldnames(tmps);
            for ii = 1:length(thesefieldnames)
                f = thesefieldnames{ii};
                if isprop(obj,f)
                    obj.(f) = tmps.(f);
                end
            end
            
            % update grid accordingly
            updateLsfInnerGrid(obj);
            
        end
        
        function saveProfile(obj,src,event)
            %save current values to a profile of the class (or a class) by name
            
            % make sure up-to-date
            updateLsf(obj);
            
            % check if values are valid?
            if ~isValid(obj)
                return
            end
            
            [filename,pathname] = uiputfile('*.mat','Select a file to save');
            if filename==0
                return
            end
            file_and_path = [pathname, filesep,filename];
            
            thesefieldnames = obj.PROPERTIES_TO_SAVE;
            
            for ii = 1:length(thesefieldnames)                
                %save back to the file:
                %this is a little tricky. We want to save the information in s as
                %variables with the name obj.PROPERTIES_TO_SAVE ('lsf' for example). So we
                %make a field of a structure with the name obj.varname and the value
                %s. Then the save command we use with the '-struct' option, which
                %takes all the fields of a struct and saves them to variables
                %with the same names in a mat file. The append option means we
                %will add new variables to the mat file and replace existing ones
                %but leave other variables unchanged.
                varname = thesefieldnames{ii};
                values = obj.(varname);
                saveStruct.(varname) = values;
                if exist(file_and_path,'file')==2
                    save(file_and_path,'-struct','saveStruct',varname,'-append');
                else
                    save(file_and_path,'-struct','saveStruct',varname);
                end

            end
            
            
        end
        
        function out = isValid(obj)
            warning('SGRLAB:NotYetImplemented','validation function not yet written. It should highlight the wrong cells and let the user fix them.');
            out = true;
        end
end
    
    
    methods (Access = public) %make private?
        function obj = createComponents(obj)
            %obj.Parent.Title = 'lsf panel';
            
            %             obj.DropDown = uidropdown(obj.Parent);
            %             obj.DropDown.Position = [86 410 100 22];
            %             obj.DropDown.Items = {'Select One'};
            %             obj.DropDown.DropDownOpeningFcn = createCallbackFcn(obj.Parent.Parent, @DropDownOpening, true);
            %             obj.DropDown.ValueChangedFcn = createCallbackFcn(obj.Parent.Parent, @DropDownValueChanged, true);
            
            %            obj.DropDownLabel = uilabel(obj.Parent);
            obj.OuterGrid = uigridlayout(obj.Parent);
            obj.OuterGrid.RowHeight = {22,22,22,'1x',22}; % 5 rows
            obj.OuterGrid.ColumnWidth = {150,'1x'}; % 2 cols
            nrows = length(obj.OuterGrid.RowHeight);
            ncols = length(obj.OuterGrid.ColumnWidth);
            
            
            %             obj.DropDownLabel = uilabel(obj.OuterGrid);
            %             obj.DropDownLabel.HorizontalAlignment = 'right';
            %             obj.DropDownLabel.Position = [0 170 66 22];
            %             obj.DropDownLabel.Text = 'FFCF';
            
            %             obj.DropDown = uidropdown(obj.Parent,...
            %             'Position',[66 170 100 22],...
            %             'Items',{'Select One'},...
            %             'DropDownOpeningFcn', @(src,event) DropDownOpening(obj,src),...
            %             'ValueChangedFcn', @(src,event) DropDownValueChanged(obj,src));
            
            obj.LoadProfileButton = uibutton(obj.OuterGrid);
            obj.LoadProfileButton.Text = 'Load profile';
            obj.LoadProfileButton.Layout.Row = 1;
            obj.LoadProfileButton.Layout.Column = 1;
            obj.LoadProfileButton.ButtonPushedFcn = @(src,event) loadProfile(obj,src,event);
            
            obj.SaveProfileButton = uibutton(obj.OuterGrid);
            obj.SaveProfileButton.Text = 'Save profile';
            obj.SaveProfileButton.Layout.Row = 2;
            obj.SaveProfileButton.Layout.Column = 1;
            obj.SaveProfileButton.ButtonPushedFcn = @(src,event) saveProfile(obj,src,event);
            
            obj.DoneButton = uibutton(obj.OuterGrid);
            obj.DoneButton.Text = 'Done';
            obj.DoneButton.Layout.Row = nrows; %last row
            obj.DoneButton.Layout.Column = 1;
            
            obj.DropDown = uidropdown(obj.OuterGrid,...
                'Items',{'Select One'},...
                'DropDownOpeningFcn', @(src,event) DropDownOpening(obj,src),...
                'ValueChangedFcn', @(src,event) DropDownValueChanged(obj,src));
            obj.DropDown.Layout.Row = 3;
            obj.DropDown.Layout.Column = 1;
            
            obj.InnerGrid = uigridlayout(obj.OuterGrid);
            obj.InnerGrid.RowHeight = {22,22};
            obj.InnerGrid.ColumnWidth = {'1x','1x'};
            obj.InnerGrid.Layout.Row = [1 nrows];
            obj.InnerGrid.Layout.Column = 2;
        end
        
        
        function delete(obj)
            delete(obj.DropDown)
            delete(obj.DropDownLabel)
        end
    end
end

%uigetfile
%uitable.ColumnEditable=true
%g = uigridlayout([size_row size_col])
%ef = uieditfield(g,'numeric');
%ef.Layout.Row = 3;
%ef.Layout.Column = 2;

% scroll function on figure...
