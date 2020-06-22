function varargout = fileviewer(varargin)
% FILEVIEWER MATLAB code for fileviewer.fig
%      FILEVIEWER, by itself, creates a new FILEVIEWER or raises the existing
%      singleton*.
%
%      H = FILEVIEWER returns the handle to a new FILEVIEWER or the handle to
%      the existing singleton*.
%
%      FILEVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FILEVIEWER.M with the given input arguments.
%
%      FILEVIEWER('Property','Value',...) creates a new FILEVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fileviewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fileviewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fileviewer

% Last Modified by GUIDE v2.5 11-Dec-2015 09:03:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fileviewer_OpeningFcn, ...
                   'gui_OutputFcn',  @fileviewer_OutputFcn, ...
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


% --- Executes just before fileviewer is made visible.
function fileviewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fileviewer (see VARARGIN)

% Choose default command line output for fileviewer
addpath ('\\128.112.86.75\lithium\IMAGING DATA\Image Analysis');
handles.output = hObject;
%files
handles.pathname = '';
handles.filename='';
handles.fullpath = '';
handles.files_in_dir='';
%image metadata
handles.stamp='';
handles.keys='';
handles.vals='';
%image
handles.Atoms = 1;
handles.Beam = 2;
handles.Dark = 3;
handles.images = '';
handles.od = '';
handles.odProcessed = '';
%display
handles.x1Disp = '';
handles.x2Disp = '';
handles.y1Disp = '';
handles.y2Disp = '';
handles.MaxOD = 1;
handles.MinOD = 0;
%last displayed state
handles.x1Prev = '';
handles.x2Prev = '';
handles.y1Prev = '';
handles.y2Prev = '';
%processing
handles.binSize = 1;
handles.RemovedLines = 0;
handles.NormOD = 0;

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes fileviewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = fileviewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
 
% --- Executes on button press in chooseImgFile.
function chooseImgFile_Callback(hObject, eventdata, handles)
% hObject    handle to chooseImgFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Next update: remove other imaging processing buttons and have all the work
%done by the update button. Add text boxes for the other guys.
start_path=handles.pathname;
[filename,pathname]=uigetfile({'*.aia';'*.fits';'All Files (*.*)'},'Choose file to display',start_path);
%write data to handles
handles.pathname=pathname;
handles.filename=filename;
file_list=[dir(fullfile(handles.pathname,'*.aia'));dir(fullfile(handles.pathname,'*.fits'))];
files_in_dir=cell(1,length(file_list));
for ii=1:length(file_list)
    files_in_dir{ii}=file_list(ii).name;
end
handles.files_in_dir=files_in_dir;
handles.fullpath=fullfile(pathname,filename);

set(handles.text2,'string',handles.pathname);
set(handles.locationText,'string',handles.filename);


%get image data
[y1Disp,y2Disp,x1Disp,x2Disp,y1Norm,y2Norm,x1Norm,x2Norm,MinOD,MaxOD,Atoms,Beam,Dark] = getImageParams(handles);
handles.y1Disp = y1Disp;
handles.y2Disp = y2Disp;
handles.x1Disp = x1Disp;
handles.x2Disp = x2Disp;
handles.y1Last = y1Disp;
handles.y2Last = y2Disp;
handles.x1Last = x1Disp;
handles.x2Last = x2Disp;
handles.y1Norm = y1Norm;
handles.y2Norm = y2Norm;
handles.x1Norm = x1Norm;
handles.x2Norm = x2Norm;
handles.MinOD = MinOD;
handles.MaxOD = MaxOD;
handles.Atoms = Atoms;
handles.Beam = Beam;
handles.Dark = Dark;
%read image file
[Images,~,stamp,vals,keys]=readimg(handles.fullpath);
%set handles
handles.images = Images;
handles.od = getOD(Images(:,:,[handles.Atoms,handles.Beam,handles.Dark]));
handles.stamp=stamp;
handles.keys=keys;
handles.vals=vals;
set(handles.uitable1,'data',transpose(cat(1,keys,num2cell(vals))));
guidata(hObject, handles);

handles.odProcessed = handles.od(handles.y1Disp:handles.y2Disp,handles.x1Disp:handles.x2Disp);
guidata(hObject, handles);

imagesc(handles.odProcessed,[handles.MinOD,handles.MaxOD])
colorbar;





function x1dispEdit_Callback(hObject, eventdata, handles)
% hObject    handle to x1dispEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x1dispEdit as text
%        str2double(get(hObject,'String')) returns contents of x1dispEdit as a double


% --- Executes during object creation, after setting all properties.
function x1dispEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x1dispEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x2dispEdit_Callback(hObject, eventdata, handles)
% hObject    handle to x2dispEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x2dispEdit as text
%        str2double(get(hObject,'String')) returns contents of x2dispEdit as a double


% --- Executes during object creation, after setting all properties.
function x2dispEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x2dispEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y1dispEdit_Callback(hObject, eventdata, handles)
% hObject    handle to y1dispEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y1dispEdit as text
%        str2double(get(hObject,'String')) returns contents of y1dispEdit as a double


% --- Executes during object creation, after setting all properties.
function y1dispEdit_CreateFcn(hObject, ~, handles)
% hObject    handle to y1dispEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y2dispEdit_Callback(hObject, eventdata, handles)
% hObject    handle to y2dispEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y2dispEdit as text
%        str2double(get(hObject,'String')) returns contents of y2dispEdit as a double


% --- Executes during object creation, after setting all properties.
function y2dispEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y2dispEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in updatebutton.
function updatebutton_Callback(hObject, eventdata, handles)
% hObject    handle to updatebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    %Figure out how to get pic=get() to work
    
    
   [y1Disp,y2Disp,x1Disp,x2Disp,y1Norm,y2Norm,x1Norm,x2Norm,MinOD,MaxOD,Atoms,Beam,Dark] = getImageParams(handles);
   handles.y1Disp = y1Disp;
   handles.y2Disp = y2Disp;
   handles.x1Disp = x1Disp;
   handles.x2Disp = x2Disp;
   handles.y1Norm = y1Norm;
   handles.y2Norm = y2Norm;
   handles.x1Norm = x1Norm;
   handles.x2Norm = x2Norm;
   handles.MinOD = MinOD;
   handles.MaxOD = MaxOD;
   handles.Atoms = Atoms;
   handles.Beam = Beam;
   handles.Dark = Dark;
   guidata(hObject, handles);
   
   
   
   %if necessary, process image
   if handles.y1Disp~=handles.y1Last || handles.y2Disp~=handles.y2Last || handles.x1Disp~=handles.x1Last || handles.x2Disp~=handles.x2Last
       
       %normalize OD
       if handles.NormOD == 1
           ImagesODCorrection = getNormOD(handles.images(:,:,[handles.Atoms,handles.Beam,handles.Dark]),handles.images(handles.y1Norm:handles.y2Norm,handles.x1Norm:handles.x2Norm,[handles.Atoms,handles.Beam,handles.Dark]));
       else
           ImagesODCorrection = handles.images(:,:,[handles.Atoms,handles.Beam,handles.Dark]); 
       end
       
       %remove stripes and crop
       if handles.RemovedLines == 1
           ImagesCropped = ImagesODCorrection(handles.y1Disp:handles.y2Disp,handles.x1Disp:handles.x2Disp,:);
           ImagesStripeCorrection = normalizeStripes(ImagesCropped,ImagesODCorrection(handles.y1Norm:handles.y2Norm,handles.x1Norm:handles.x2Norm,:));
       else
           ImagesStripeCorrection = ImagesODCorrection(handles.y1Disp:handles.y2Disp,handles.x1Disp:handles.x2Disp,:);
       end
       
       odProc = getOD(ImagesStripeCorrection);
       
       if handles.binSize~=1
           handles.odProcessed = binImg(odProc,handles.binSize,handles.binSize);
       else
           handles.odProcessed = odProc;
       end
       
       guidata(hObject, handles);
   else
   end
   
   
   imagesc(handles.odProcessed,[handles.MinOD,handles.MaxOD])
   colorbar;


   % --- Executes on button press in previousImg.
   function previousImg_Callback(hObject, eventdata, handles)
       % hObject    handle to previousImg (see GCBO)
       % eventdata  reserved - to be defined in a future version of MATLAB
       % handles    structure with handles and user data (see GUIDATA)
       [~,index]=ismember(handles.filename,handles.files_in_dir);
       if index>1
           filename=handles.files_in_dir{index-1};
       elseif index==1
           filename=handles.files_in_dir{length(handles.files_in_dir)};
       else
       end
       
       handles.fullpath=fullfile(handles.pathname,filename);
       handles.filename=filename;
       
       set(handles.text2,'string',handles.pathname);
       set(handles.locationText,'string',handles.filename);
       
       
       %get image data
       [y1Disp,y2Disp,x1Disp,x2Disp,y1Norm,y2Norm,x1Norm,x2Norm,MinOD,MaxOD,Atoms,Beam,Dark] = getImageParams(handles);
       handles.y1Disp = y1Disp;
       handles.y2Disp = y2Disp;
       handles.x1Disp = x1Disp;
       handles.x2Disp = x2Disp;
       handles.y1Last = y1Disp;
       handles.y2Last = y2Disp;
       handles.x1Last = x1Disp;
       handles.x2Last = x2Disp;
       handles.y1Norm = y1Norm;
       handles.y2Norm = y2Norm;
       handles.x1Norm = x1Norm;
       handles.x2Norm = x2Norm;
       handles.MinOD = MinOD;
       handles.MaxOD = MaxOD;
       handles.Atoms = Atoms;
       handles.Beam = Beam;
       handles.Dark = Dark;
       [Images,~,stamp,vals,keys]=readimg(handles.fullpath);
       %set handles
       handles.images = Images;
       handles.od = getOD(Images(:,:,[handles.Atoms,handles.Beam,handles.Dark]));
       handles.stamp=stamp;
       handles.keys=keys;
       handles.vals=vals;
       set(handles.uitable1,'data',transpose(cat(1,keys,num2cell(vals))));
       
       handles.odProcessed = handles.od(handles.y1Disp:handles.y2Disp,handles.x1Disp:handles.x2Disp);
       guidata(hObject, handles);
       
       imagesc(handles.odProcessed,[handles.MinOD,handles.MaxOD])
       colorbar;


       % --- Executes on button press in nextImg.
       function nextImg_Callback(hObject, eventdata, handles)
           % hObject    handle to nextImg (see GCBO)
           % eventdata  reserved - to be defined in a future version of MATLAB
           % handles    structure with handles and user data (see GUIDATA)
           [~,index]=ismember(handles.filename,handles.files_in_dir);
           if index<length(handles.files_in_dir)
               filename=handles.files_in_dir{index+1};
           elseif index==length(handles.files_in_dir)
               filename=handles.files_in_dir{1};
           else
           end
           
           handles.fullpath=fullfile(handles.pathname,filename);
           handles.filename=filename;
           
           
           set(handles.text2,'string',handles.pathname);
           set(handles.locationText,'string',handles.filename);
           
           
           %get image data
           %get image data
           [y1Disp,y2Disp,x1Disp,x2Disp,y1Norm,y2Norm,x1Norm,x2Norm,MinOD,MaxOD,Atoms,Beam,Dark] = getImageParams(handles);
           handles.y1Disp = y1Disp;
           handles.y2Disp = y2Disp;
           handles.x1Disp = x1Disp;
           handles.x2Disp = x2Disp;
           handles.y1Last = y1Disp;
           handles.y2Last = y2Disp;
           handles.x1Last = x1Disp;
           handles.x2Last = x2Disp;
           handles.y1Norm = y1Norm;
           handles.y2Norm = y2Norm;
           handles.x1Norm = x1Norm;
           handles.x2Norm = x2Norm;
           handles.MinOD = MinOD;
           handles.MaxOD = MaxOD;
           handles.Atoms = Atoms;
           handles.Beam = Beam;
           handles.Dark = Dark;
           [Images,~,stamp,vals,keys]=readimg(handles.fullpath);
           %set handles
           handles.images = Images;
           handles.od = getOD(Images(:,:,[handles.Atoms,handles.Beam,handles.Dark]));
           handles.stamp=stamp;
           handles.keys=keys;
           handles.vals=vals;
           set(handles.uitable1,'data',transpose(cat(1,keys,num2cell(vals))));
           guidata(hObject, handles);
           handles.odProcessed = handles.od(handles.y1Disp:handles.y2Disp,handles.x1Disp:handles.x2Disp);
           guidata(hObject, handles);
           
           imagesc(handles.odProcessed,[handles.MinOD,handles.MaxOD])
           colorbar;



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double



% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function binEdit_Callback(hObject, eventdata, handles)
% hObject    handle to binEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of binEdit as text
%        str2double(get(hObject,'String')) returns contents of binEdit as a double



% --- Executes during object creation, after setting all properties.
function binEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to binEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in removeLinesButton.
    function pushbutton5_Callback(hObject, eventdata, handles)
        % hObject    handle to removeLinesButton (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)
        
        handles.RemovedLines = 1;
        
        [y1Disp,y2Disp,x1Disp,x2Disp,y1Norm,y2Norm,x1Norm,x2Norm,MinOD,MaxOD,Atoms,Beam,Dark] = getImageParams(handles);
        handles.y1Disp = y1Disp;
        handles.y2Disp = y2Disp;
        handles.x1Disp = x1Disp;
        handles.x2Disp = x2Disp;
        handles.y1Norm = y1Norm;
        handles.y2Norm = y2Norm;
        handles.x1Norm = x1Norm;
        handles.x2Norm = x2Norm;
        handles.MinOD = MinOD;
        handles.MaxOD = MaxOD;
        handles.Atoms = Atoms;
        handles.Beam = Beam;
        handles.Dark = Dark;
        guidata(hObject, handles);
        
        %normalize OD
        if handles.NormOD == 1
            ImagesODCorrection = getNormOD(handles.images(:,:,[handles.Atoms,handles.Beam,handles.Dark]),handles.images(handles.y1Norm:handles.y2Norm,handles.x1Norm:handles.x2Norm,[handles.Atoms,handles.Beam,handles.Dark]));
        else
            ImagesODCorrection = handles.images(:,:,[handles.Atoms,handles.Beam,handles.Dark]);
        end
        
        %remove stripes and crop
        
        ImagesCropped = ImagesODCorrection(handles.y1Disp:handles.y2Disp,handles.x1Disp:handles.x2Disp,:);
        ImagesStripeCorrection = normalizeStripes(ImagesCropped,ImagesODCorrection(handles.y1Norm:handles.y2Norm,handles.x1Norm:handles.x2Norm,:));
        
        odProc = getOD(ImagesStripeCorrection);
        
        if handles.binSize~=1
            handles.odProcessed = binImg(odProc,handles.binSize,handles.binSize);
        else
            handles.odProcessed = odProc;
        end
        
        guidata(hObject, handles);
        
        
        imagesc(handles.odProcessed,[handles.MinOD,handles.MaxOD])
        colorbar;



% --- Executes on button press in removeLinesButton.
function removeLinesButton_Callback(hObject, eventdata, handles)
% hObject    handle to removeLinesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double



% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double



% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function odMinEdit_Callback(hObject, eventdata, handles)
% hObject    handle to odMinEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of odMinEdit as text
%        str2double(get(hObject,'String')) returns contents of odMinEdit as a double


% --- Executes during object creation, after setting all properties.
function odMinEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to odMinEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function odMaxEdit_Callback(hObject, eventdata, handles)
% hObject    handle to odMaxEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of odMaxEdit as text
%        str2double(get(hObject,'String')) returns contents of odMaxEdit as a double


% --- Executes during object creation, after setting all properties.
function odMaxEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to odMaxEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in normODbutton.
function normODbutton_Callback(hObject, eventdata, handles)
% hObject    handle to normODbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.NormOD = 1;
[y1Disp,y2Disp,x1Disp,x2Disp,y1Norm,y2Norm,x1Norm,x2Norm,MinOD,MaxOD,Atoms,Beam,Dark] = getImageParams(handles);
handles.y1Disp = y1Disp;
handles.y2Disp = y2Disp;
handles.x1Disp = x1Disp;
handles.x2Disp = x2Disp;
handles.y1Norm = y1Norm;
handles.y2Norm = y2Norm;
handles.x1Norm = x1Norm;
handles.x2Norm = x2Norm;
handles.MinOD = MinOD;
handles.MaxOD = MaxOD;
handles.Atoms = Atoms;
handles.Beam = Beam;
handles.Dark = Dark;
guidata(hObject, handles);
%normalize OD
ImagesODCorrection = getNormOD(handles.images(:,:,[handles.Atoms,handles.Beam,handles.Dark]),handles.images(handles.y1Norm:handles.y2Norm,handles.x1Norm:handles.x2Norm,[handles.Atoms,handles.Beam,handles.Dark]));

%remove stripes and crop
if handles.RemovedLines == 1
    ImagesCropped = ImagesODCorrection(handles.y1Disp:handles.y2Disp,handles.x1Disp:handles.x2Disp,:);
    ImagesStripeCorrection = normalizeStripes(ImagesCropped,ImagesODCorrection(handles.y1Norm:handles.y2Norm,handles.x1Norm:handles.x2Norm,:));
else
    ImagesStripeCorrection = ImagesODCorrection(handles.y1Disp:handles.y2Disp,handles.x1Disp:handles.x2Disp,:);
end

odProc = getOD(ImagesStripeCorrection);

if handles.binSize~=1
    handles.odProcessed = binImg(odProc,handles.binSize,handles.binSize);
else
    handles.odProcessed = odProc;
end
imagesc(handles.odProcessed,[handles.MinOD,handles.MaxOD])
colorbar;



function x1normEdit_Callback(hObject, eventdata, handles)
% hObject    handle to x1normEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x1normEdit as text
%        str2double(get(hObject,'String')) returns contents of x1normEdit as a double


% --- Executes during object creation, after setting all properties.
function x1normEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x1normEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x2normEdit_Callback(hObject, eventdata, handles)
% hObject    handle to x2normEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x2normEdit as text
%        str2double(get(hObject,'String')) returns contents of x2normEdit as a double


% --- Executes during object creation, after setting all properties.
function x2normEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x2normEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y1normEdit_Callback(hObject, eventdata, handles)
% hObject    handle to y1normEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y1normEdit as text
%        str2double(get(hObject,'String')) returns contents of y1normEdit as a double


% --- Executes during object creation, after setting all properties.
function y1normEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y1normEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y2normEdit_Callback(hObject, eventdata, handles)
% hObject    handle to y2normEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y2normEdit as text
%        str2double(get(hObject,'String')) returns contents of y2normEdit as a double


% --- Executes during object creation, after setting all properties.
function y2normEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y2normEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function atomsEdit_Callback(hObject, eventdata, handles)
% hObject    handle to atomsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of atomsEdit as text
%        str2double(get(hObject,'String')) returns contents of atomsEdit as a double


% --- Executes during object creation, after setting all properties.
function atomsEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to atomsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function beamEdit_Callback(hObject, eventdata, handles)
% hObject    handle to beamEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beamEdit as text
%        str2double(get(hObject,'String')) returns contents of beamEdit as a double


% --- Executes during object creation, after setting all properties.
function beamEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beamEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function darkEdit_Callback(hObject, eventdata, handles)
% hObject    handle to darkEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of darkEdit as text
%        str2double(get(hObject,'String')) returns contents of darkEdit as a double


% --- Executes during object creation, after setting all properties.
function darkEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to darkEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 

function [y1Disp,y2Disp,x1Disp,x2Disp,y1Norm,y2Norm,x1Norm,x2Norm,MinOD,MaxOD,Atoms,Beam,Dark] = getImageParams(handles)
        %function should be called before plotting. Updates all parameters
        %set by text boxes.
        
        %Image region
        y1strDisp=get(handles.y1dispEdit, 'String');
        y2strDisp=get(handles.y2dispEdit, 'String');
        x1strDisp=get(handles.x1dispEdit, 'String');
        x2strDisp=get(handles.x2dispEdit, 'String');
        
        if (~isempty(y1strDisp)&&~isempty(y2strDisp)&&~isempty(x1strDisp)&&~isempty(x2strDisp))
            y1Disp=str2double(y1strDisp);
            y2Disp=str2double(y2strDisp);
            x1Disp=str2double(x1strDisp);
            x2Disp=str2double(x2strDisp);
        else
            y1Disp=1;
            y2Disp=size(handles.od,1);
            x1Disp=1;
            x2Disp=size(handles.od,2);
        end
        
        %normalization region
        y1strNorm=get(handles.y1normEdit, 'String');
        y2strNorm=get(handles.y2normEdit, 'String');
        x1strNorm=get(handles.x1normEdit, 'String');
        x2strNorm=get(handles.x2normEdit, 'String');
        
        if (~isempty(y1strNorm)&&~isempty(y2strNorm)&&~isempty(x1strNorm)&&~isempty(x2strNorm))
            y1Norm=str2double(y1strNorm);
            y2Norm=str2double(y2strNorm);
            x1Norm=str2double(x1strNorm);
            x2Norm=str2double(x2strNorm);
        else
            y1Norm=1;
            y2Norm=size(handles.od,1);
            x1Norm=1;
            x2Norm=size(handles.od,2);
        end
        
        
        %od bounds
        MinODStr = get(handles.odMinEdit,'String');
        MaxODStr = get(handles.odMaxEdit,'String');
        if (~isempty(MinODStr)&&~isempty(MaxODStr))
            MinOD = str2double(MinODStr);
            MaxOD = str2double(MaxODStr);
        else
            MinOD = -0.2;
            MaxOD = 1;
        end
        
        AtomsStr = get(handles.atomsEdit,'String');
        BeamStr = get(handles.beamEdit,'String');
        DarkStr = get(handles.darkEdit,'String');
        if (~isempty(AtomsStr)&&~isempty(BeamStr)&&~isempty(DarkStr))
            Atoms = str2double(AtomsStr);
            Beam = str2double(BeamStr);
            Dark = str2double(DarkStr);
        else
            Atoms = 1;
            Beam = 2;
            Dark = 3;
        end


% --- Executes on button press in binbutton.
function binbutton_Callback(hObject, eventdata, handles)
% hObject    handle to binbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   binSizeStr = get(handles.binEdit,'String');
   if ~isempty(binSizeStr)
       handles.binSize = str2double(binSizeStr);
   else
       handles.binSize = 1;
   end
   guidata(hObject, handles)
   
   [y1Disp,y2Disp,x1Disp,x2Disp,y1Norm,y2Norm,x1Norm,x2Norm,MinOD,MaxOD,~,~,~] = getImageParams(handles);
   handles.y1Disp = y1Disp;
   handles.y2Disp = y2Disp;
   handles.x1Disp = x1Disp;
   handles.x2Disp = x2Disp;
   handles.y1Norm = y1Norm;
   handles.y2Norm = y2Norm;
   handles.x1Norm = x1Norm;
   handles.x2Norm = x2Norm;
   handles.MinOD = MinOD;
   handles.MaxOD = MaxOD;
      
   handles.odProcessed = binImg(handles.odProcessed,handles.binSize,handles.binSize);
   guidata(hObject, handles);
   
   imagesc(handles.odProcessed,[handles.MinOD,handles.MaxOD])
   colorbar;
   
    
   
