function varargout = ManualSpikeCurateGUI(varargin)
% MANUALSPIKECURATEGUI MATLAB code for ManualSpikeCurateGUI.fig
%      MANUALSPIKECURATEGUI, by itself, creates a new MANUALSPIKECURATEGUI or raises the existing
%      singleton*.
%
%      H = MANUALSPIKECURATEGUI returns the handle to a new MANUALSPIKECURATEGUI or the handle to
%      the existing singleton*.
%
%      MANUALSPIKECURATEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MANUALSPIKECURATEGUI.M with the given input arguments.
%
%      MANUALSPIKECURATEGUI('Property','Value',...) creates a new MANUALSPIKECURATEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ManualSpikeCurateGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ManualSpikeCurateGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ManualSpikeCurateGUI

% Last Modified by GUIDE v2.5 21-Oct-2023 15:16:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ManualSpikeCurateGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ManualSpikeCurateGUI_OutputFcn, ...
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


% --- Executes just before ManualSpikeCurateGUI is made visible.
function ManualSpikeCurateGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ManualSpikeCurateGUI (see VARARGIN)

% Choose default command line output for ManualSpikeCurateGUI
handles.output = hObject;
initialData = NaN(600,1);
set(handles.uitable1,'data',initialData)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ManualSpikeCurateGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ManualSpikeCurateGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uigetfile('*.mat','Select Calcium .mat File');
load([path file])
handles.PSTHhitspks = Spikes.PSTH.hit.spks;
handles.PSTHhitFR = Spikes.PSTH.hit.spkRates;
handles.PSTHmissspks = Spikes.PSTH.miss.spks;
handles.PSTHmissFR = Spikes.PSTH.miss.spkRates;
handles.ROI = length(handles.PSTHhitspks);
% kernal for spike ploting
handles.kernal = [0 1 1 0;0 1 1 0;0 1 1 0];

%set the slider value to 1 and set the max/min values depending on numbre
%of ROIs
sliderMin = 1;
sliderMax = handles.ROI;
set(handles.slider1,'SliderStep',[1, 1] / (sliderMax - sliderMin),'max',handles.ROI,'min',1,'Value',1)
x=get(handles.slider1,'Value');

%set the current component string textbox to 1
set(handles.edit1,'string',num2str(x));

% show rasterPlot
spkTemp = conv2(handles.PSTHhitspks{x},handles.kernal,'valid');
spkTemp2 = conv2(handles.PSTHmissspks{x},handles.kernal,'valid');
imagesc(handles.axes1,spkTemp);colormap(flip(gray)),caxis([0 1]);
imagesc(handles.axes3,spkTemp2);colormap(flip(gray)),caxis([0 1]);
set(handles.axes1,'YTickLabel','','XTickLabel','','YTick','','XTick','')
% show rasterPlot
bar(handles.axes2,handles.PSTHhitFR(x,:));
bar(handles.axes4,handles.PSTHmissFR(x,:));
%if the variable badComponents exists, update these Bad Components table
if(exist('badComponents','var'))
    initialData = NaN(600,1);
    initialData(1:length(badComponents)) = badComponents;
    set(handles.uitable1,'data',initialData); 
end

%update handles
guidata(hObject, handles);


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
x=get(handles.slider1,'Value');
x = ceil(x);
%update the text box to new slider value
set(handles.edit1,'string',num2str(x));

%replot spikes and FR
spkTemp = conv2(handles.PSTHhitspks{x},handles.kernal,'valid');
spkTemp2 = conv2(handles.PSTHmissspks{x},handles.kernal,'valid');
imagesc(handles.axes1,spkTemp);colormap(flip(gray)),caxis([0 1]);
imagesc(handles.axes3,spkTemp2);colormap(flip(gray)),caxis([0 1]);
set(handles.axes1,'YTickLabel','','XTickLabel','','YTick','','XTick','')
% show rasterPlot
bar(handles.axes2,handles.PSTHhitFR(x,:));
bar(handles.axes4,handles.PSTHmissFR(x,:));


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Mark as bad component
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get the current component selection
currentComponent = str2num(get(handles.edit1,'string'));
%get the current data in the table
tableData=get(handles.uitable1,'data');
%find all NaN
nans = find(isnan(tableData));
%set the first NaN to the current component selection
tableData(nans(1)) = currentComponent;
%reset the table data
set(handles.uitable1,'data',tableData);
% Increment slider
% set(handles.slider1,'Value',currentComponent+1);
guidata(hObject, handles);



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
%get the input value
x=str2num(get(handles.edit1,'string'));
set(handles.slider1,'Value',x);
spkTemp = conv2(handles.PSTHhitspks{x},handles.kernal,'valid');
spkTemp2 = conv2(handles.PSTHmissspks{x},handles.kernal,'valid');
imagesc(handles.axes1,spkTemp);colormap(flip(gray)),caxis([0 1]);
imagesc(handles.axes3,spkTemp2);colormap(flip(gray)),caxis([0 1]);
set(handles.axes1,'YTickLabel','','XTickLabel','','YTick','','XTick','')
% show rasterPlot
bar(handles.axes2,handles.PSTHhitFR(x,:));
bar(handles.axes4,handles.PSTHmissFR(x,:));
% Plots spike and FR


guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% -- Send to workspace
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%get the table data with all bad components
tableData=get(handles.uitable1,'data');
%find all NaN
nans = find(isnan(tableData));
%set all NaNs to empty
tableData(nans) = [];
%send table data to workspace
assignin('base','goodSpkComponents',tableData)
