function varargout = Helen_Miranda_GUI_Input(varargin)
% Helen_Miranda_GUI_Input MATLAB code for Helen_Miranda_GUI_Input.fig
%      Helen_Miranda_GUI_Input, by itself, creates a new Helen_Miranda_GUI_Input or raises the existing
%      singleton*.
%
%      H = Helen_Miranda_GUI_Input returns the handle to a new Helen_Miranda_GUI_Input or the handle to
%      the existing singleton*.
%
%      Helen_Miranda_GUI_Input('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in Helen_Miranda_GUI_Input.M with the given input arguments.
%
%      Helen_Miranda_GUI_Input('Property','Value',...) creates a new Helen_Miranda_GUI_Input or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the Helen_Miranda_GUI_Input before Helen_Miranda_GUI_Input_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Helen_Miranda_GUI_Input_OpeningFcn via varargin.
%
%      *See Helen_Miranda_GUI_Input Options on GUIDE's Tools menu.  Choose "Helen_Miranda_GUI_Input allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Helen_Miranda_GUI_Input

% Last Modified by GUIDE v2.5 25-Jul-2023 20:45:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Helen_Miranda_GUI_Input_OpeningFcn, ...
                   'gui_OutputFcn',  @Helen_Miranda_GUI_Input_OutputFcn, ...
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


% --- Executes just before Helen_Miranda_GUI_Input is made visible.
function Helen_Miranda_GUI_Input_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Helen_Miranda_GUI_Input (see VARARGIN)

% Choose default command line output for Helen_Miranda_GUI_Input
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Helen_Miranda_GUI_Input wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Helen_Miranda_GUI_Input_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --- Executes during object creation, after setting all properties.
function Thickness_of_Plate_TextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Thickness_of_Plate_TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function Elastic_Modulus_TextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Elastic_Modulus_TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RUN_Button.
function RUN_Button_Callback(hObject, eventdata, handles)
% hObject    handle to RUN_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

thickness_of_plate = evalin('base','thickness_of_plate');
nodof = evalin('base','nodof');
dim = evalin('base','dim');
Boundary_Conditions = evalin('base','Boundary_Conditions');
External_Load = evalin('base','External_Load');
nodal_coordinate_values = evalin('base','nodal_coordinate_values');
nodal_connectivity_values = evalin('base','nodal_connectivity_values');
Length_of_Element = evalin('base','Length_of_Element');
Width_of_element = evalin('base','Width_of_element');
Element_Type = evalin('base','Element_Type');

analysis_main = Analysis(Length_of_Element,Width_of_element ...
    ,thickness_of_plate,nodof,dim,Boundary_Conditions ...
    ,External_Load,nodal_coordinate_values,nodal_connectivity_values,Element_Type);

analysis_main.Elastic_Modulus = evalin('base','Elastic_Modulus');
analysis_main.poissons_ratio = evalin('base','Possions_Ratio');

analysis_main.sigma_t = evalin('base','sigma_t');
analysis_main.Percentage_Limit_Tension = evalin('base','Percentage_Limit_Tension');
analysis_main.Percentage_Limit_Shear = evalin('base','Percentage_Limit_Shear');
analysis_main.tau = evalin('base','tau');
analysis_main.Percentage_Fracture_Energy_Traction = evalin('base','Percentage_Fracture_Energy_Traction');
analysis_main.d_u = evalin('base','d_u');
analysis_main.Percentage_Ultimate_Deformation_Maximum_Step = evalin('base','Percentage_Ultimate_Deformation_Maximum_Step');
analysis_main.tol = evalin('base','tol');
analysis_main.gamma = evalin('base','weigth_density');

analysis_main.is_dx_enabled = evalin('base','is_dx_enabled');
analysis_main.is_dy_enabled = evalin('base','is_dy_enabled');
analysis_main.is_dz_enabled = evalin('base','is_dz_enabled');
analysis_main.is_rx_enabled = evalin('base','is_rx_enabled');
analysis_main.is_ry_enabled = evalin('base','is_ry_enabled');
analysis_main.is_rz_enabled = evalin('base','is_rz_enabled');

analysis_main.Engine(analysis_main);

% Helen_Miranda_4_Quadilateral_No_Inputs


% --- Executes during object creation, after setting all properties.
function Possions_Ratio_TextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Possions_Ratio_TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function Length_of_Element_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Length_of_Element_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function width_of_element_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to width_of_element_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function NXE_Textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NXE_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function NXY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NXY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function nodof_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nodof_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CoordinatesBrowseButton.
function CoordinatesBrowseButton_Callback(hObject, eventdata, handles)

            [FileName, Path] = uigetfile('*.txt'); %Will display File Explorer to Choose a text file
            Complete_path = strcat(Path,'\',FileName); %Concatenates string to use in declaring variable to workspace
            nodal_coordinate_values = load(Complete_path); %Declaring Variable to workspace
            assignin('base','nodal_coordinate_values',nodal_coordinate_values);
            
            set(handles.CoordinatesFileName_Label,'String',FileName);
            
            clear Complete_path FileName Path

% hObject    handle to CoordinatesBrowseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function CoordinatesFileName_Label_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BoundaryConditionFileName_Label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in NodalConnectivityBrowseButton.
function NodalConnectivityBrowseButton_Callback(hObject, eventdata, handles)
% hObject    handle to NodalConnectivityBrowseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

            [FileName, Path] = uigetfile('*.txt'); %Will display File Explorer to Choose a text file
            Complete_path = strcat(Path,'\',FileName); %Concatenates string to use in declaring variable to workspace
            nodal_connectivity_values = load(Complete_path); %Declaring Variable to workspace
            assignin('base','nodal_connectivity_values',nodal_connectivity_values);
            
            set(handles.NodalConnectivityFileName_Label,'String',FileName);
            
            clear Complete_path FileName Path



% --- Executes during object creation, after setting all properties.
function NodalConnectivityFileName_Label_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NodalConnectivityFileName_Label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BoundaryConditionBrowseButton.
function BoundaryConditionBrowseButton_Callback(hObject, eventdata, handles)
% hObject    handle to BoundaryConditionBrowseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

            [FileName, Path] = uigetfile('*.txt'); %Will display File Explorer to Choose a text file
            Complete_path = strcat(Path,'\',FileName); %Concatenates string to use in declaring variable to workspace
            Boundry_Conditions = load(Complete_path); %Declaring Variable to workspace
            assignin('base','Boundary_Conditions',Boundry_Conditions);
            set(handles.BoundaryConditionFileName_Label,'String',FileName);
            
            clear Complete_path FileName Path



% --- Executes during object creation, after setting all properties.
function BoundaryConditionFileName_Label_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BoundaryConditionFileName_Label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in ExternalLoadBrowseButton.
function ExternalLoadBrowseButton_Callback(hObject, eventdata, handles)
% hObject    handle to ExternalLoadBrowseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
            [FileName, Path] = uigetfile('*.txt'); %Will display File Explorer to Choose a text file
            Complete_path = strcat(Path,'\',FileName); %Concatenates string to use in declaring variable to workspace
            External_Load = load(Complete_path); %Declaring Variable to workspace
            assignin('base','External_Load',External_Load);
            set(handles.ExternalLoadFileNameLabel,'String',FileName);
            
            
            clear Complete_path FileName Path


% --- Executes during object creation, after setting all properties.
function ExternalLoadFileNameLabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ExternalLoadFileNameLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SaveandPlotMesh.
function SaveandPlotMesh_Callback(hObject, eventdata, handles)
% hObject    handle to SaveandPlotMesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Global Inputs
thickness_of_plate = str2double(get(handles.Thickness_of_Plate_TextBox,'String'));
assignin('base','thickness_of_plate',thickness_of_plate);

Elastic_Modulus = str2double(get(handles.Elastic_Modulus_TextBox,'String'));
assignin('base','Elastic_Modulus',Elastic_Modulus);

Possions_Ratio = str2double(get(handles.Possions_Ratio_TextBox,'String'));
assignin('base','Possions_Ratio',Possions_Ratio);


% Element Inputs
Length_of_Element = str2double(get(handles.Length_of_Element_Textbox,'String'));
assignin('base','Length_of_Element',Length_of_Element);

Width_of_element = str2double(get(handles.width_of_element_textbox,'String'));
assignin('base','Width_of_element',Width_of_element);

NXE = str2double(get(handles.NXE_Textbox,'String'));
assignin('base','NXE',NXE);

NYE = str2double(get(handles.NXY,'String'));
assignin('base','NYE',NYE);

nodof = str2double(get(handles.nodof_textbox,'String'));
assignin('base','nodof',nodof);

Element_Type = str2double(get(handles.Element_Type_Edit,'String'));
assignin('base','Element_Type',Element_Type);

ngpb = str2double(get(handles.ngpb,'String'));
assignin('base','ngpb',ngpb);

ngps = str2double(get(handles.ngps,'String'));
assignin('base','ngps',ngps);

weigth_density = str2double(get(handles.Weight_Density_TextBox,'String'));
assignin('base','weigth_density',weigth_density);

if get(handles.dim_checkbox,'Value')
    dim = 2;
    assignin('base','dim',dim);
else
    dim = 3;
    assignin('base','dim',dim);
end

sigma_t = str2double(get(handles.sigma_t,'String'));
assignin('base','sigma_t',sigma_t);

Percentage_Limit_Tension = str2double(get(handles.Percentage_Limit_Tension,'String'));
assignin('base','Percentage_Limit_Tension',Percentage_Limit_Tension);

tau = str2double(get(handles.tau,'String'));
assignin('base','tau',tau);

Percentage_Limit_Shear = str2double(get(handles.Percentage_Limit_Shear,'String'));
assignin('base','Percentage_Limit_Shear',Percentage_Limit_Shear);

Gf_t = str2double(get(handles.Gf_t,'String'));
assignin('base','Gf_t',Gf_t);

Percentage_Fracture_Energy_Traction = str2double(get(handles.Percentage_Fracture_Energy_Traction,'String'));
assignin('base','Percentage_Fracture_Energy_Traction',Percentage_Fracture_Energy_Traction);

d_u = str2double(get(handles.d_u,'String'));
assignin('base','d_u',d_u);

Percentage_Ultimate_Deformation_Maximum_Step  = str2double(get(handles.Percentage_Ultimate_Deformation_Maximum_Step ,'String'));
assignin('base','Percentage_Ultimate_Deformation_Maximum_Step',Percentage_Ultimate_Deformation_Maximum_Step);

tol  = str2double(get(handles.tol ,'String'));
assignin('base','tol',tol );


dx_enabled_bool = true;
dy_enabled_bool = true;
dz_enabled_bool = true;
rx_enabled_bool = true;
rx_enabled_bool = true;
rx_enabled_bool = true;

%checking if the user wants to include dx in the analysis : Dx
if handles.dx_enabled.Value
    dx_enabled_bool = true;
    assignin('base','is_dx_enabled',dx_enabled_bool);
else
    dx_enabled_bool = false;
    assignin('base','is_dx_enabled',dx_enabled_bool);
end


%checking if the user wants to include dy in the analysis : Dy
if handles.dy_enabled.Value
    dy_enabled_bool = true;
    assignin('base','is_dy_enabled',dy_enabled_bool);
else
    dy_enabled_bool = false;
    assignin('base','is_dy_enabled',dy_enabled_bool);
end

%checking if the user wants to include dz in the analysis: Dz
if handles.dz_enabled.Value
    dz_enabled_bool = true;
    assignin('base','is_dz_enabled',dz_enabled_bool);
else
    dz_enabled_bool = false;
    assignin('base','is_dz_enabled',dz_enabled_bool);
end

%checking if the user wants to include rx in the analysis: Rx
if handles.rx_enabled.Value
    rx_enabled_bool = true;
    assignin('base','is_rx_enabled',rx_enabled_bool);
else
    rx_enabled_bool = false;
    assignin('base','is_rx_enabled',rx_enabled_bool);
end

%checking if the user wants to include rx in the analysis: Ry
if handles.ry_enabled.Value
    ry_enabled_bool = true;
    assignin('base','is_ry_enabled',ry_enabled_bool);
else
    ry_enabled_bool = false;
    assignin('base','is_ry_enabled',ry_enabled_bool);
end

%checking if the user wants to include rx in the analysis: Rz
if handles.rz_enabled.Value
    rz_enabled_bool = true;
    assignin('base','is_rz_enabled',rz_enabled_bool);
else
    rz_enabled_bool = false;
    assignin('base','is_rz_enabled',rz_enabled_bool);
end











% --- Executes during object deletion, before destroying properties.
function CoordinatesFileName_Label_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to CoordinatesFileName_Label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in dim_checkbox.
function dim_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to dim_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% Hint: get(hObject,'Value') returns toggle state of dim_checkbox


% --- Executes during object creation, after setting all properties.
function dim_checkbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dim_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Making the Mesh Plot

nodal_coordinate_values = evalin('base', 'nodal_coordinate_values');
nodal_connectivity_values = evalin('base', 'nodal_connectivity_values');

PlotMesh(nodal_coordinate_values,nodal_connectivity_values);



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Element_Type_Edit.
function Element_Type_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Element_Type_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Element_Type_Edit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Element_Type_Edit


% --- Executes during object creation, after setting all properties.
function Element_Type_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Element_Type_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ngpb_Callback(hObject, eventdata, handles)
% hObject    handle to ngpb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ngpb as text
%        str2double(get(hObject,'String')) returns contents of ngpb as a double


% --- Executes during object creation, after setting all properties.
function ngpb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ngpb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ngps_Callback(hObject, eventdata, handles)
% hObject    handle to ngps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ngps as text
%        str2double(get(hObject,'String')) returns contents of ngps as a double


% --- Executes during object creation, after setting all properties.
function ngps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ngps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dx_enabled.
function dx_enabled_Callback(hObject, eventdata, handles)
% hObject    handle to dx_enabled (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dx_enabled


% --- Executes during object creation, after setting all properties.
function dx_enabled_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dx_enabled (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in dy_enabled.
function dy_enabled_Callback(hObject, eventdata, handles)
% hObject    handle to dy_enabled (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dy_enabled


% --- Executes on button press in dz_enabled.
function dz_enabled_Callback(hObject, eventdata, handles)
% hObject    handle to dz_enabled (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dz_enabled


% --- Executes on button press in rx_enabled.
function rx_enabled_Callback(hObject, eventdata, handles)
% hObject    handle to rx_enabled (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rx_enabled


% --- Executes on button press in ry_enabled.
function ry_enabled_Callback(hObject, eventdata, handles)
% hObject    handle to ry_enabled (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ry_enabled


% --- Executes on button press in rz_enabled.
function rz_enabled_Callback(hObject, eventdata, handles)
% hObject    handle to rz_enabled (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rz_enabled



function Weight_Density_Label_Callback(hObject, eventdata, handles)
% hObject    handle to Weight_Density_Label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Weight_Density_Label as text
%        str2double(get(hObject,'String')) returns contents of Weight_Density_Label as a double


% --- Executes during object creation, after setting all properties.
function Weight_Density_Label_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Weight_Density_Label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Weight_Density_TextBox_Callback(hObject, eventdata, handles)
% hObject    handle to Weight_Density_TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Weight_Density_TextBox as text
%        str2double(get(hObject,'String')) returns contents of Weight_Density_TextBox as a double


% --- Executes during object creation, after setting all properties.
function Weight_Density_TextBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Weight_Density_TextBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigma_t_Callback(hObject, eventdata, handles)
% hObject    handle to sigma_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigma_t as text
%        str2double(get(hObject,'String')) returns contents of sigma_t as a double


% --- Executes during object creation, after setting all properties.
function sigma_t_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigma_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Percentage_Limit_Tension_Callback(hObject, eventdata, handles)
% hObject    handle to Percentage_Limit_Tension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Percentage_Limit_Tension as text
%        str2double(get(hObject,'String')) returns contents of Percentage_Limit_Tension as a double


% --- Executes during object creation, after setting all properties.
function Percentage_Limit_Tension_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Percentage_Limit_Tension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit27 as text
%        str2double(get(hObject,'String')) returns contents of edit27 as a double


% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tau_Callback(hObject, eventdata, handles)
% hObject    handle to tau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tau as text
%        str2double(get(hObject,'String')) returns contents of tau as a double


% --- Executes during object creation, after setting all properties.
function tau_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit29_Callback(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit29 as text
%        str2double(get(hObject,'String')) returns contents of edit29 as a double


% --- Executes during object creation, after setting all properties.
function edit29_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Percentage_Limit_Shear_Callback(hObject, eventdata, handles)
% hObject    handle to Percentage_Limit_Shear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Percentage_Limit_Shear as text
%        str2double(get(hObject,'String')) returns contents of Percentage_Limit_Shear as a double


% --- Executes during object creation, after setting all properties.
function Percentage_Limit_Shear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Percentage_Limit_Shear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit31_Callback(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit31 as text
%        str2double(get(hObject,'String')) returns contents of edit31 as a double


% --- Executes during object creation, after setting all properties.
function edit31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Gf_t_Callback(hObject, eventdata, handles)
% hObject    handle to Gf_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gf_t as text
%        str2double(get(hObject,'String')) returns contents of Gf_t as a double


% --- Executes during object creation, after setting all properties.
function Gf_t_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gf_t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit33_Callback(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit33 as text
%        str2double(get(hObject,'String')) returns contents of edit33 as a double


% --- Executes during object creation, after setting all properties.
function edit33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Percentage_Fracture_Energy_Traction_Callback(hObject, eventdata, handles)
% hObject    handle to Percentage_Fracture_Energy_Traction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Percentage_Fracture_Energy_Traction as text
%        str2double(get(hObject,'String')) returns contents of Percentage_Fracture_Energy_Traction as a double


% --- Executes during object creation, after setting all properties.
function Percentage_Fracture_Energy_Traction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Percentage_Fracture_Energy_Traction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit35_Callback(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit35 as text
%        str2double(get(hObject,'String')) returns contents of edit35 as a double


% --- Executes during object creation, after setting all properties.
function edit35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function d_u_Callback(hObject, eventdata, handles)
% hObject    handle to d_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d_u as text
%        str2double(get(hObject,'String')) returns contents of d_u as a double


% --- Executes during object creation, after setting all properties.
function d_u_CreateFcn(hObject, eventdata, handles)
% hObject    handle to d_u (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit37_Callback(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit37 as text
%        str2double(get(hObject,'String')) returns contents of edit37 as a double


% --- Executes during object creation, after setting all properties.
function edit37_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Percentage_Ultimate_Deformation_Maximum_Step_Callback(hObject, eventdata, handles)
% hObject    handle to Percentage_Ultimate_Deformation_Maximum_Step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Percentage_Ultimate_Deformation_Maximum_Step as text
%        str2double(get(hObject,'String')) returns contents of Percentage_Ultimate_Deformation_Maximum_Step as a double


% --- Executes during object creation, after setting all properties.
function Percentage_Ultimate_Deformation_Maximum_Step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Percentage_Ultimate_Deformation_Maximum_Step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit39_Callback(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit39 as text
%        str2double(get(hObject,'String')) returns contents of edit39 as a double


% --- Executes during object creation, after setting all properties.
function edit39_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tol_Callback(hObject, eventdata, handles)
% hObject    handle to tol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tol as text
%        str2double(get(hObject,'String')) returns contents of tol as a double


% --- Executes during object creation, after setting all properties.
function tol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Length_of_Element_Textbox_Callback(hObject, eventdata, handles)
% hObject    handle to Length_of_Element_Textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of Length_of_Element_Textbox as text
%        str2double(get(hObject,'String')) returns contents of Length_of_Element_Textbox as a double
