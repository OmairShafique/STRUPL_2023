function varargout = Call_Me(varargin)
% CALL_ME MATLAB code for Call_Me.fig
%      CALL_ME, by itself, creates a new CALL_ME or raises the existing
%      singleton*.
%
%      H = CALL_ME returns the handle to a new CALL_ME or the handle to
%      the existing singleton*.
%
%      CALL_ME('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALL_ME.M with the given input arguments.
%
%      CALL_ME('Property','Value',...) creates a new CALL_ME or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Call_Me_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Call_Me_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Call_Me

% Last Modified by GUIDE v2.5 25-Jul-2023 21:03:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Call_Me_OpeningFcn, ...
                   'gui_OutputFcn',  @Call_Me_OutputFcn, ...
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


% --- Executes just before Call_Me is made visible.
function Call_Me_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Call_Me (see VARARGIN)

% Choose default command line output for Call_Me
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Call_Me wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Call_Me_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Call_This_Button.
function Call_This_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Call_This_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
