function varargout = load_2DIR(varargin)
% LOAD_2DIR M-file for load_2DIR.fig
%      LOAD_2DIR, by itself, creates a new LOAD_2DIR or raises the existing
%      singleton*.
%
%      H = LOAD_2DIR returns the handle to a new LOAD_2DIR or the handle to
%      the existing singleton*.
%
%      LOAD_2DIR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOAD_2DIR.M with the given input arguments.
%
%      LOAD_2DIR('Property','Value',...) creates a new LOAD_2DIR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before load_2DIR_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to load_2DIR_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help load_2DIR

% Last Modified by GUIDE v2.5 04-Jun-2011 23:55:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @load_2DIR_OpeningFcn, ...
                   'gui_OutputFcn',  @load_2DIR_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before load_2DIR is made visible.
function load_2DIR_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to load_2DIR (see VARARGIN)

% Choose default command line output for load_2DIR
handles.output = hObject;

load_dir_list(handles);
handles.raw = construct2d();
handles.raw.zeropad = 1;
handles.raw.undersampling = 1;
handles.raw.fft_type = 'sgrsfft';
handles.raw.freq_units = 'cm-1';
handles.raw.time_units = 'fs';

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes load_2DIR wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = load_2DIR_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in dir_list.
function dir_list_Callback(hObject, eventdata, handles)
% hObject    handle to dir_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns dir_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dir_list

index_selected = get(hObject,'Value');
handles.dir_names = get(hObject,'String');
dir_sel = handles.dir_names{index_selected};
handles.meta = is2DIR(dir_sel, handles);
if handles.meta.vers;
    handles.raw.meta_text = handles.meta.text;
    handles.raw.basename = char(dir_sel);
    set(handles.meta_text, 'String',handles.meta.text);
    set(handles.var_name_edit, 'String',handles.raw.basename);
else
    if strcmp(get(handles.figure1,'SelectionType'),'open');
        cd(dir_sel);
        set(hObject, 'Value', 1);
        set(handles.meta_text,'String','Not a data directory.');
        load_dir_list(handles);
    end
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function dir_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dir_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Populate directory list
function load_dir_list(handles)
% load_dir_list populates the list box Directory List
% with all the directories according to the chosen menu.
% ie, all data directories, photon counter data directories, 
% or spectral data directories. Other directories or files are omitted
% excepted for '..' to go up one or '.' to refresh.

total_dir = dir;
sorted_dir = [];

for i = 1:size(total_dir,1);
    if total_dir(i).isdir
        sorted_dir = [sorted_dir;total_dir(i)];
    end
end
[names, index] = sortrows({sorted_dir.name}');
handles.dir_names = names;
handles.dir_index = index;
guidata(handles.figure1,handles);
set(handles.dir_list,'String',handles.dir_names,'Value',1);
set(handles.figure1,'Name',pwd);
guidata(handles.figure1,handles);

% --- Determines if dir has 2DIR data
function out = is2DIR(temp_dir,handles)
% this program is a bit of a cludge. It only checks to see if one of the 
% files inside the directory is named "*_meta.txt" It does not verify the
% data any further.
% temp_dir - directory string
% out - structure containing struct with .vers, .text, .lines


dir_name = [char(temp_dir),'/*_meta.txt'];
if size(dir(dir_name),1);
    try meta_path = dir([char(temp_dir),'/*_meta.txt']);
        meta_file = fopen([char(temp_dir),'/',meta_path.name]);
        meta = textscan(meta_file, '%s', 'whitespace', '');
        fclose(meta_file);
        out.text = char(meta{1});
        lines = textscan(out.text, '%s','whitespace','\n');
        out.lines = lines{1};
        words = textscan(out.text, '%s','whitespace');
        out.vers = str2num(char(words{1}(2)));
        
    catch
        out.vers = 0;
    end
else
    out.vers = 0;
end



% --- Parse Meta Data
function out = meta_parser(handles)
% frewind(handles.meta_file);
% meta = textscan(handles.meta_file, '%s', 'whitespace', '\n');
% meta = meta{1};
% line = textscan(meta{1}, '%s');
% handles.meta.ver = str2num(char(line{1}(2)));

if handles.meta.vers == 1.1;
    temp_line = key_val('Steps', handles.meta.lines);
    handles.raw.Steps = str2num(char(temp_line{1}(2)));
    temp_line = key_val('StepSize', handles.meta.lines);
    handles.raw.StepSize = str2num(char(temp_line{1}(2)));
    temp_line = key_val('Shots', handles.meta.lines);
    handles.raw.Shots = str2num(char(temp_line{1}(2)));
    temp_line = key_val('T2', handles.meta.lines);
    handles.raw.t2 = str2num(char(temp_line{1}(2)));
    temp_line = key_val('T3', handles.meta.lines);
    handles.raw.t3 = str2num(char(temp_line{1}(2)));
    temp_line = key_val('Phase', handles.meta.lines);
    degree = char(176); %% cludge
    handles.raw.phase = sscanf(char(temp_line{1}(2)),...
        ['%d' degree]);
    temp_line = key_val('Center wavelength', handles.meta.lines);
    handles.raw.centerwavelength = str2num(char(temp_line{1}(3)));
    temp_line = key_val('Center wavenumber', handles.meta.lines);
    handles.raw.centerwavenumber = str2num(char(temp_line{1}(3)));
else
    % other code goes here for other versions
end
out = handles;

function var_name_edit_Callback(hObject, eventdata, handles)
% hObject    handle to var_name_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of var_name_edit as text
%        str2double(get(hObject,'String')) returns contents of var_name_edit as a double


% --- Executes during object creation, after setting all properties.
function var_name_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to var_name_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zero_pad_edit_Callback(hObject, eventdata, handles)
% hObject    handle to zero_pad_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zero_pad_edit as text
%        str2double(get(hObject,'String')) returns contents of zero_pad_edit as a double


% --- Executes during object creation, after setting all properties.
function zero_pad_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zero_pad_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function under_sample_edit_Callback(hObject, eventdata, handles)
% hObject    handle to under_sample_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of under_sample_edit as text
%        str2double(get(hObject,'String')) returns contents of under_sample_edit as a double


% --- Executes during object creation, after setting all properties.
function under_sample_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to under_sample_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in fft_type_menu.
function fft_type_menu_Callback(hObject, eventdata, handles)
% hObject    handle to fft_type_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fft_type_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        fft_type_menu



% --- Executes during object creation, after setting all properties.
function fft_type_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fft_type_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_var_button.
function save_var_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_var_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
       
try
    handles.var_name = get(handles.var_name_edit,'String');
    assignin('base',handles.var_name,[]); % Save structure to workspace
catch exception
    errordlg(sprintf('"%s" is not valid.\nTry a different variable name.',...
        handles.var_name),...
        'Bad input', 'modal');
    rethrow(exception);
end
handles = make_var(handles);
assignin('base',handles.var_name,handles.raw);

function raw = make_var(handles)
file_dirs = get(handles.dir_list, 'String');
dir_val = get(handles.dir_list, 'Value');
dir_sel = file_dirs{dir_val};
handles.raw.zeropad = str2num(get(handles.zero_pad_edit,'String'));
handles.raw.undersampling = str2num(get(handles.under_sample_edit,'String'));
fft_type_list = {'sgrsfft','fft'};
handles.raw.fft_type = fft_type_list{get(handles.fft_type_menu,'Value')};
if handles.meta.vers;
    handles = meta_parser(handles);
    handles.raw.basename = handles.var_name;
    handles.raw.data.R = load([char(dir_sel),'/',char(dir_sel),'.dat']);
    handles.raw.data.NR = load([char(dir_sel),'/',char(dir_sel),'_NR.dat']);
    handles.raw.freq = handles.raw.data.R(1,2:end);
    handles.raw.freq_units = 'cm-1';
    handles.raw.time = handles.raw.data.R(2:end,1);
    handles.raw.time_units = 'fs';
    handles.raw.R1 = handles.raw.data.R(2:end,2:end);
    handles.raw.R2 = handles.raw.data.NR(2:end,2:end);
    handles.raw = freq2d(handles.raw);
    handles.raw.w3 = handles.raw.freq - 0;
%     raw = absorptive2d(raw,...
%         'phase', handles.meta.Phase,...
%         'zeropad',(handles.meta.zeropad * length(raw.time)),...
%         'range', [raw.freq(1),raw.freq(end)],...
%         'fft_type',handles.meta.fft_type);
    raw = handles;
    
end



function val = key_val(key, meta_data)
val = 0;
idx = strfind(meta_data, key);
for i = 1:size(meta_data,1);
    if isempty(idx{i})
    else
        meta_data{i};
        val = textscan(meta_data{i},'%s');
        break
    end
end
