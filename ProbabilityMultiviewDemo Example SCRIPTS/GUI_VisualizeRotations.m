function varargout = GUI_VisualizeRotations(varargin)
% GUI_VISUALIZEROTATIONS visualize 3x3 matrices and check for valid
% rotations.
%	GUI_VISUALIZEROTATIONS allows users to input 3x3 matrices and
%	visualize them in 3D space using two relative coordinate frames.
%	Indicators also show users is the matrix is a valid rotation matrix
%	to help show the associated properties.
% 
%   NOTES:
%     (1) This GUI should be run from the MATLAB command window. Do not use
%         the "Run" button to run this code.
%     (2) This GUI requires the Transformation Toolbox that is included in
%         the installation of the ScorBot Toolbox. Please goto 
%         https://www.usna.edu/Users/weapsys/kutzer/_Code-Development/ScorBot_Toolbox.php
%         and follow the download and install instructions prior to running
%         this GUI.
%
%   See also: triad isSO 
%
%   M. Kutzer, 13Aug2015, USNA

% This code was initialized using the MATLAB GUIDE function

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_VisualizeRotations_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_VisualizeRotations_OutputFcn, ...
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


% --- Executes just before GUI_VisualizeRotations is made visible.
function GUI_VisualizeRotations_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_VisualizeRotations (see VARARGIN)

% TODO - Eliminate global variables
global h1 body vertices

% Check for ScorBot Toolbox
switch exist('triad','file')
    case 2
        % At least one file from the toolbox is installed!
    otherwise
        fprintf(2,'Please install the ScorBot toolbox before running this code.\n\n');
        web('https://www.usna.edu/Users/weapsys/kutzer/_Code-Development/ScorBot_Toolbox.php','-browser');
end

% Choose default command line output for GUI_VisualizeRotations
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Convert to 3D view and turn grid on
view(handles.axs,3);
grid(handles.axs,'on');

% Label axes
xlabel(handles.axs,'x_0');
ylabel(handles.axs,'y_0');
zlabel(handles.axs,'z_0');

% Create frame 0
%h0 = triad('parent',handles.axs,'axislabels',{'x_0','y_0','z_0'},'linewidth',1);
h0 = triad('parent',handles.axs,'linewidth',1,'LineStyle','--');
% Create frame 1
h1 = triad('parent',handles.axs,'axislabels',{'x_1','y_1','z_1'},'linewidth',1.5);

% Load Wall-E
open('Wall-E.fig');
ff = gcf;
set(ff,'Visible','off');
aa = get(ff,'Children');
body = get(aa,'Children');
set(body,'Parent',handles.axs,'FaceAlpha',0.7);
close(ff);
addSingleLight(handles.axs);
vertices = body.Vertices;

rotate3d(handles.axs);

updateFrame(handles);
% UIWAIT makes GUI_VisualizeRotations wait for user response (see UIRESUME)
% uiwait(handles.GUI_VisualizeRotations);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_VisualizeRotations_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles;


function R_11_Callback(hObject, eventdata, handles)
% hObject    handle to R_11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R_11 as text
%        str2double(get(hObject,'String')) returns contents of R_11 as a double
updateFrame(handles);


% --- Executes during object creation, after setting all properties.
function R_11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R_11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function R_12_Callback(hObject, eventdata, handles)
% hObject    handle to R_12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R_12 as text
%        str2double(get(hObject,'String')) returns contents of R_12 as a double
updateFrame(handles);


% --- Executes during object creation, after setting all properties.
function R_12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R_12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function R_13_Callback(hObject, eventdata, handles)
% hObject    handle to R_13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R_13 as text
%        str2double(get(hObject,'String')) returns contents of R_13 as a double
updateFrame(handles);


% --- Executes during object creation, after setting all properties.
function R_13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R_13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function R_21_Callback(hObject, eventdata, handles)
% hObject    handle to R_21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R_21 as text
%        str2double(get(hObject,'String')) returns contents of R_21 as a double
updateFrame(handles);


% --- Executes during object creation, after setting all properties.
function R_21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R_21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function R_22_Callback(hObject, eventdata, handles)
% hObject    handle to R_22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R_22 as text
%        str2double(get(hObject,'String')) returns contents of R_22 as a double
updateFrame(handles);


% --- Executes during object creation, after setting all properties.
function R_22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R_22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function R_23_Callback(hObject, eventdata, handles)
% hObject    handle to R_23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R_23 as text
%        str2double(get(hObject,'String')) returns contents of R_23 as a double
updateFrame(handles);


% --- Executes during object creation, after setting all properties.
function R_23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R_23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function R_31_Callback(hObject, eventdata, handles)
% hObject    handle to R_31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R_31 as text
%        str2double(get(hObject,'String')) returns contents of R_31 as a double
updateFrame(handles);


% --- Executes during object creation, after setting all properties.
function R_31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R_31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function R_32_Callback(hObject, eventdata, handles)
% hObject    handle to R_32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R_32 as text
%        str2double(get(hObject,'String')) returns contents of R_32 as a double
updateFrame(handles);


% --- Executes during object creation, after setting all properties.
function R_32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R_32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function R_33_Callback(hObject, eventdata, handles)
% hObject    handle to R_33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of R_33 as text
%        str2double(get(hObject,'String')) returns contents of R_33 as a double
updateFrame(handles);


% --- Executes during object creation, after setting all properties.
function R_33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to R_33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in OrthogonalColumns.
function OrthogonalColumns_Callback(hObject, eventdata, handles)
% hObject    handle to OrthogonalColumns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OrthogonalColumns
flag = get(handles.OrthogonalColumns,'Value');
switch flag
    case 1
        R = getRotation(handles);
        for j = 1:3
            R_mag(j) = norm(R(:,j));
            R_hat(:,j) = R(:,j)./R_mag(j);
        end
        
        R_star{1} = diag([ 1, 1, 1]);
        R_star{2} = diag([-1, 1, 1]);
        %R_star{3} = diag([ 1,-1, 1]);
        %R_star{4} = diag([ 1, 1,-1]);
        %R_star{5} = diag([ 1,-1,-1]);
        %R_star{6} = diag([-1, 1,-1]);
        %R_star{7} = diag([ 1,-1,-1]);
        %R_star{8} = diag([-1,-1,-1]);
        
        for i = 1:numel(R_star)
            p = [R_hat,zeros(3,1)];
            q = [R_star{i},zeros(3,1)];
            [H_q2p{i},err(i)] = pointsToSE3(q,p);
        end
        [~,idx] = min(err);
        R = H_q2p{idx}(1:3,1:3)*diag(R_mag)*R_star{idx};
        
        setRotation(handles,R);
        updateFrame(handles);
    case 0
        set(handles.OrthogonalColumns,'Value',1);
end


% --- Executes on button press in OrthogonalRows.
function OrthogonalRows_Callback(hObject, eventdata, handles)
% hObject    handle to OrthogonalRows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OrthogonalRows
flag = get(handles.OrthogonalRows,'Value');
switch flag
    case 1
        R = getRotation(handles);
        for i = 1:3
            R_mag(i) = norm(R(i,:));
            R_hat(i,:) = R(i,:)./R_mag(i);
        end
        
        R_star{1} = diag([ 1, 1, 1]);
        R_star{2} = diag([-1, 1, 1]);
        %R_star{3} = diag([ 1,-1, 1]);
        %R_star{4} = diag([ 1, 1,-1]);
        %R_star{5} = diag([ 1,-1,-1]);
        %R_star{6} = diag([-1, 1,-1]);
        %R_star{7} = diag([ 1,-1,-1]);
        %R_star{8} = diag([-1,-1,-1]);
        
        for i = 1:numel(R_star)
            p = [transpose(R_hat),zeros(3,1)];
            q = [R_star{i},zeros(3,1)];
            [H_q2p{i},err(i)] = pointsToSE3(q,p);
        end
        [~,idx] = min(err);
        R = transpose(H_q2p{idx}(1:3,1:3)*diag(R_mag)*R_star{idx});
        
        setRotation(handles,R);
        updateFrame(handles);
    case 0
        set(handles.OrthogonalRows,'Value',1);
end

% --- Executes on button press in UnitVectorRows.
function UnitVectorRows_Callback(hObject, eventdata, handles)
% hObject    handle to UnitVectorRows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UnitVectorRows
flag = get(handles.UnitVectorRows,'Value');
switch flag
    case 1
        R = getRotation(handles);
        for i = 1:3
            R(i,:) = R(i,:)./norm(R(i,:));
        end
        setRotation(handles,R);
        updateFrame(handles);
    case 0
        set(handles.UnitVectorRows,'Value',1);
end

% --- Executes on button press in UnitVectorColumns.
function UnitVectorColumns_Callback(hObject, eventdata, handles)
% hObject    handle to UnitVectorColumns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of UnitVectorColumns
flag = get(handles.UnitVectorColumns,'Value');
switch flag
    case 1
        R = getRotation(handles);
        for j = 1:3
            R(:,j) = R(:,j)./norm(R(:,j));
        end
        setRotation(handles,R);
        updateFrame(handles);
    case 0
        set(handles.UnitVectorColumns,'Value',1);
end

% --- Executes on button press in Determinant.
function Determinant_Callback(hObject, eventdata, handles)
% hObject    handle to Determinant (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Determinant
flag = get(handles.Determinant,'Value');
switch flag
    case 1
        R = getRotation(handles);
        d = det(R);
        s = sign(d) * ( abs(d)^(1/3) );
        R = R*diag(repmat(1/s,3,1));
        setRotation(handles,R);
        updateFrame(handles);
    case 0
        set(handles.Determinant,'Value',1);
end

% --- Executes on button press in InverseTranspose.
function InverseTranspose_Callback(hObject, eventdata, handles)
% hObject    handle to InverseTranspose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of InverseTranspose
flag = get(handles.InverseTranspose,'Value');
switch flag
    case 1
        R = getRotation(handles);
        for j = 1:3
            R_mag(j) = norm(R(:,j));
            R_hat(:,j) = R(:,j)./R_mag(j);
        end
        
        R_star{1} = diag([ 1, 1, 1]);
        R_star{2} = diag([-1, 1, 1]);
        %R_star{3} = diag([ 1,-1, 1]);
        %R_star{4} = diag([ 1, 1,-1]);
        %R_star{5} = diag([ 1,-1,-1]);
        %R_star{6} = diag([-1, 1,-1]);
        %R_star{7} = diag([ 1,-1,-1]);
        %R_star{8} = diag([-1,-1,-1]);
        
        for i = 1:numel(R_star)
            p = [R_hat,zeros(3,1)];
            q = [R_star{i},zeros(3,1)];
            [H_q2p{i},err(i)] = pointsToSE3(q,p);
        end
        [~,idx] = min(err);
        R = H_q2p{idx}(1:3,1:3)*sign(diag(R_mag))*R_star{idx};
        
        setRotation(handles,R);
        updateFrame(handles);
    case 0
        set(handles.InverseTranspose,'Value',1);
end


% --- Executes on button press in ValidRotation.
function ValidRotation_Callback(hObject, eventdata, handles)
% hObject    handle to ValidRotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ValidRotation
flag = get(handles.ValidRotation,'Value');
switch flag
    case 1
        R = getRotation(handles);
        for j = 1:3
            R_mag(j) = norm(R(:,j));
            R_hat(:,j) = R(:,j)./R_mag(j);
        end
        
        R_star{1} = diag([ 1, 1, 1]);
        R_star{2} = diag([-1, 1, 1]);
        %R_star{3} = diag([ 1,-1, 1]);
        %R_star{4} = diag([ 1, 1,-1]);
        %R_star{5} = diag([ 1,-1,-1]);
        %R_star{6} = diag([-1, 1,-1]);
        %R_star{7} = diag([ 1,-1,-1]);
        %R_star{8} = diag([-1,-1,-1]);
        
        for i = 1:numel(R_star)
            p = [R_hat,zeros(3,1)];
            q = [R_star{i},zeros(3,1)];
            [H_q2p{i},err(i)] = pointsToSE3(q,p);
        end
        [~,idx] = min(err);
        R = H_q2p{idx}(1:3,1:3);
        
        setRotation(handles,R);
        updateFrame(handles);
    case 0
        set(handles.ValidRotation,'Value',1);
end


% --- Internal function to update frame and check matrix
function updateFrame(handles)

% TODO - Eliminate global variables
global h1 ZERO body vertices

% Update rotation matrix
R = getRotation(handles);
setRotation(handles, R);

% Get handles for h1
kids = get(h1,'Children');
% Get x,y,z axis handles
x =  findobj(kids,'tag','X-Axis');
y =  findobj(kids,'tag','Y-Axis');
z =  findobj(kids,'tag','Z-Axis');
% Update x,y,z axis data
set(x,'xData',[0,R(1,1)],'yData',[0,R(2,1)],'zData',[0,R(3,1)]);
set(y,'xData',[0,R(1,2)],'yData',[0,R(2,2)],'zData',[0,R(3,2)]);
set(z,'xData',[0,R(1,3)],'yData',[0,R(2,3)],'zData',[0,R(3,3)]);
% Get x,y,z axis labels
xlbl =  findobj(kids,'tag','X-Label');
ylbl =  findobj(kids,'tag','Y-Label');
zlbl =  findobj(kids,'tag','Z-Label');
% Update x,y,z axis label positions
set(xlbl,'Position',R(:,1)');
set(ylbl,'Position',R(:,2)');
set(zlbl,'Position',R(:,3)');

% Update STL
try
    v = R*transpose(vertices);
    body.Vertices = transpose(v);
    set(body,'Visible','on');
catch
    set(body,'Visible',off);
end

% Check rotation matrix
ZERO = 0.005; %2e3*eps(class(R));

% (1) Orthogonal columns
flag(1) = true;
% Check for zero length vectors (they cannot be orthogonal)
for j = 1:3
    if norm(R(:,j)) <= ZERO
        flag(1) = false;
    end
end
% Check for orthogonality
if max( abs( dot(R(:,1),R(:,2)) ) ) > ZERO
    flag(1) = false;
end
if max( abs( dot(R(:,1),R(:,3)) ) ) > ZERO
    flag(1) = false;
end
if max( abs( dot(R(:,2),R(:,3)) ) ) > ZERO
    flag(1) = false;
end

if flag(1)
    set(handles.OrthogonalColumns,'Value',1);
else
    set(handles.OrthogonalColumns,'Value',0);
end

% (2) Orthogonal rows
flag(2) = true;
% Check for zero length vectors (they cannot be orthogonal)
for i = 1:3
    if norm(R(i,:)) <= ZERO
        flag(2) = false;
    end
end
% Check for orthogonality
if max( abs( dot(R(1,:),R(2,:)) ) ) > ZERO
    flag(2) = false;
end
if max( abs( dot(R(1,:),R(3,:)) ) ) > ZERO
    flag(2) = false;
end
if max( abs( dot(R(2,:),R(3,:)) ) ) > ZERO
    flag(2) = false;
end

if flag(2)
    set(handles.OrthogonalRows,'Value',1);
else
    set(handles.OrthogonalRows,'Value',0);
end

% (3) Unit vector columns
flag(3) = true;
for j = 1:3
    if abs(norm(R(:,j)) - 1) > ZERO
        flag(3) = false;
    end
end

if flag(3)
    set(handles.UnitVectorColumns,'Value',1);
else
    set(handles.UnitVectorColumns,'Value',0);
end

% (4) Unit vector rows
flag(4) = true;
for i = 1:3
    if abs(norm(R(i,:)) - 1) > ZERO
        flag(4) = false;
    end
end

if flag(4)
    set(handles.UnitVectorRows,'Value',1);
else
    set(handles.UnitVectorRows,'Value',0);
end

% (5) Determinant = 1
flag(5) = true;
if abs(det(R) - 1) > ZERO
    flag(5) = false;
end

if flag(5)
    set(handles.Determinant,'Value',1);
else
    set(handles.Determinant,'Value',0);
end

% (6) Inverse = Transpose
flag(6) = true;
% Check if matrix can be inverted
if abs( det(R) ) < ZERO
    flag(6) = false;
else
    % If it can be inverted, check the inverse and transpose
    Z = minv(R) - transpose(R);
    if max( max( abs(Z) ) ) > ZERO
        flag(6) = false;
    end
end

if flag(6)
    set(handles.InverseTranspose,'Value',1);
else
    set(handles.InverseTranspose,'Value',0);
end

% (7) Valid Rotation
if sum(flag) == numel(flag)
    set(handles.ValidRotation,'Value',1);
else
    set(handles.ValidRotation,'Value',0);
end

drawnow

% --- Internal function to get rotation
function R = getRotation(handles)

% Update rotation matrix
for i = 1:3
    for j = 1:3
        str = get(handles.( sprintf('R_%d%d',i,j) ), 'String');
        try
            eval( sprintf( 'R(%d,%d) = %s;',i,j,str) );
            %set(handles.( sprintf('R_%d%d',i,j) ), 'String', sprintf('%.3f',R(i,j)) );
        catch
            R(i,j) = NaN;
            %set(handles.( sprintf('R_%d%d',i,j) ), 'String','NaN');
        end
    end
end

% --- Internal function to set rotation
function setRotation(handles, R)

% Set rotation matrix in GUI
for i = 1:3
    for j = 1:3
        set(handles.( sprintf('R_%d%d',i,j) ), 'String', sprintf('%.3f',R(i,j)) );
    end
end