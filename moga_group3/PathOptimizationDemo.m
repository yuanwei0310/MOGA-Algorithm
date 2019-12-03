function varargout = PathOptimizationDemo(varargin)
% PATHOPTIMIZATIONDEMO MATLAB code for PathOptimizationDemo.fig
% This is an path-planning optimization demo using
% 
% Optimization Toolbox
% Global Optimization Toolbox
% 
% The objective is to find the optimal path through a randomly generated
% vector field of wind values.
% 
% The blue curve can be moved manually to try out different paths.
% The optimization toolbox can me called in 2 ways:
% 
% Locally: Start FMINCON from the local manual path.
% Globally: Use MultiStart (with 10 runs) to find a globally optimal path. This
% uses PCT if workers are available.
% 
% -------Parameters-------
% Airspeed: Speed of the airplane flying the wind field
% Points: Number of control points in the line
% Interpolation: Interpolation method of the control points (using INTERP1)
% Random Seed: Specifies the seed for the random number generator. Leave empty for a random field.
% Wind Rougness: A parameter from 0 to 1 determining how wavy to make the wind
% Field Size: Size of the vector field
%
% Copyright (c) 2012, MathWorks, Inc. 
%

% Edit the above text to modify the response to help PathOptimizationDemo

% Last Modified by GUIDE v2.5 05-Apr-2012 16:47:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PathOptimizationDemo_OpeningFcn, ...
                   'gui_OutputFcn',  @PathOptimizationDemo_OutputFcn, ...
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


% --- Executes just before PathOptimizationDemo is made visible.
function PathOptimizationDemo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PathOptimizationDemo (see VARARGIN)
% Choose default command line output for PathOptimizationDemo
set(handles.AxesMain,'fontsize',8);

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PathOptimizationDemo wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%Initialize the GUI with the New Wind Field button callback
pushbutton_NewField_Callback(hObject, eventdata, handles)


% --- Outputs from this function are returned to the command line.
function varargout = PathOptimizationDemo_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_NewField.
% Clears the data and generates a new random wind field
function pushbutton_NewField_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_NewField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateOptText(handles,[],[],[]); %Clear the "optimal" box
seed = str2double(get(handles.editSeed,'String')); %Set the random seed
if isfinite(seed)
    rng(typecast(single(seed),'uint32'));
else
    rng('shuffle');
end
delete(get(handles.AxesMain,'children'));
rng(123)
SZX = round(max(str2double(get(handles.edit_FieldX,'string')), 10));
SZY = round(max(str2double(get(handles.edit_FieldY,'string')), 10));
SZX = min(SZX,300); %Set some sort of upper bound...
SZY = min(SZY,300);
set(handles.edit_FieldX,'string',num2str(SZX));
set(handles.edit_FieldY,'string',num2str(SZY));

windFineness = get(handles.edit_WindRoughness,'string');
windFineness = str2double(windFineness);
windFineness = max(windFineness,0);
windFineness = min(windFineness,1);
set(handles.edit_WindRoughness,'string',num2str(windFineness));

%Make the wind fields
Ix = makeRandomWindField(windFineness,SZX,SZY);
Ix = Ix/(max(abs(Ix(:)))) * (50+150*rand);
Iy = makeRandomWindField(windFineness,SZX,SZY);
Iy = Iy/(max(abs(Iy(:)))) * (50+150*rand);

[Xg,Yg] = meshgrid(0:SZX,0:SZY);
%Penalty represents the wind component in the direction of the end point
L = (sqrt((Xg-SZX).^2 + (Yg-SZY/2).^2));
Penalty = (SZX-Xg)./L .* Ix +  Iy.*(SZY/2-Yg)./L; 
Penalty(~isfinite(Penalty)) = nan; 
IM = imagesc(Penalty); % This will be the background for the vector field
set(IM,'Xdata',[0 SZX]);
set(IM,'Ydata',[0 SZY]);
xlim([0 SZX]);
ylim([0 SZY]);
hold on;
colormap(interp1([0,1,2],[1 0 0; 1 1 1; 0 1 0],0:0.01:2));
%hsv = rgb2hsv(colormap);
%colormap((hsv2rgb(hsv*[1 0 0; 0 0.5 0; 0 0 1])));
%caxis(max(abs(caxis))*[-1 1]);
caxis(max(abs(Penalty(:)))*[-1 1]);

quiver(Xg,Yg,Ix,Iy,'k'); % Draw the vector field
h_colorbar = colorbar;
title(h_colorbar,'Tailwind')
xlabel(h_colorbar,'Headwind')
axis equal tight xy;
plot([0 SZX],[SZY SZY]/2,'k.','markersize',30);
drawnow;

% Store everything inside "handles"
handles.Ix = Ix;
handles.Iy = Iy;
handles.SZX = SZX;
handles.SZY = SZY;
handles.lineInteractiveCurve = [];
handles.h_ctrlpts = [];
handles.IntegrationFineness = 101; %Number of steps in the line integral

guidata(hObject, handles);

% Create the blue interactive line, and do the straight line path
makeInteractiveLine(handles);

function makeInteractiveLine(handles)
%This function sets up the interactive manual path.
npts = str2double(get(handles.edit_NUMPOINTS,'string'));
SZX = handles.SZX;
SZY = handles.SZY;
Ix = handles.Ix;
Iy = handles.Iy;

%Delete the current line if it exists
if ishandle(handles.lineInteractiveCurve)
    delete(handles.lineInteractiveCurve);
    delete(nonzeros(handles.h_ctrlpts));
end

METHOD = get(handles.popupmenu_METHOD,'value');
METHOD = getMETHODS(METHOD);

p0 = linspace(0,SZX,npts+2)'; % Linearly spaced in X
p0 = [p0 0*p0+SZY/2]; % Halfway up the figure in Y
hold all;

% Interpolate and draw the line
XYdata = PointsToCurve(p0,METHOD,handles.IntegrationFineness);
lineInteractiveCurve = plot(XYdata(:,1),XYdata(:,2),'color',[.5 .5 1],'linewidth',2);
set(lineInteractiveCurve,'hittest','off'); % We don't need this to be selectable

%Store the initial points, and the curve handle
handles.p0 = p0;
handles.lineInteractiveCurve = lineInteractiveCurve;

% These will be the draggable control points
for n = 2:(size(p0,1)-1)
    h_ctrlpts(n) = plot(p0(n,1),p0(n,2),'o','color',[.5 .5 1],'markerfacecolor',[.5 .5 1]);
end
h_ends = plot(p0([1 end],1),p0([1 end],2),'o','color',[.5 .5 1],'markerfacecolor',[.5 .5 1]);

handles.h_ctrlpts = h_ctrlpts;

guidata(handles.figure1, handles);

% Do the calculation for the straight line path
V_PLANESPEED = str2double(get(handles.edit_PLANESPEED,'string'));
[t_straightLine dist_straightLine] = ...
    calcTime(XYdata,Ix,Iy,V_PLANESPEED,'linear','initial',handles);
updateManualText(handles,dist_straightLine/t_straightLine,dist_straightLine,t_straightLine);
title(handles.AxesMain,sprintf('Straight Line Time: %d hours, %.1f minutes\n',floor(t_straightLine),rem(t_straightLine,1)*60));

% Set the control points' button down function (make them draggable)
set(nonzeros(h_ctrlpts),'buttondownfcn',{@ButtonDownOnLine,h_ctrlpts,handles});


function XYdata = PointsToCurve(p,METHOD,fineness)
% Interpolate the curve based on the control points and METHOD
% The number of points interpolated is defined in "fineness"
nP = size(p,1);
XYdata = [interp1(1:nP,p(:,1),linspace(1,nP,fineness)',METHOD,'extrap') interp1(1:nP,p(:,2),linspace(1,nP,fineness)',METHOD,'extrap')];


function ButtonDownOnLine(~,~,h_ctrlpts,handles)
% When we click on a control point, we set up the dragging action here.
% We modify the figure's windowbuttonmotionfcn.

npts = str2double(get(handles.edit_NUMPOINTS,'string'));
METHOD = get(handles.popupmenu_METHOD,'value');

METHOD = getMETHODS(METHOD);

V_PLANESPEED = str2double(get(handles.edit_PLANESPEED,'string'));
hp = h_ctrlpts(hittest == h_ctrlpts);
set(hp,'markerfacecolor','k');
set(gcbf,'windowbuttonupfcn',{@resetDrag,hp})
set(gcbf,'windowbuttonmotionfcn',{@onDrag, hp,h_ctrlpts,METHOD,V_PLANESPEED,npts,handles})

function resetDrag(~,~,hp)
% Turn of dragging when the mouse button is released
try, set(hp,'markerfacecolor',[0.5 0.5 1]); end %Sometimes this hp is not valid
set(gcbf,'windowbuttonmotionfcn','');

function onDrag(~,~,hp,h_ctrlpts, METHOD,V_PLANESPEED,npts,handles)
% This function drags the control point to the new location
% and updates the time calculation and displays it.

ax = get(hp,'parent');
cp = get(ax,'currentpoint');
% Don't let the control points go outside the box:
cp = max(cp,0);
set(hp,'xdata',min(cp(1),handles.SZX),'ydata',min(cp(3),handles.SZY));

CalculateManualValues(handles,h_ctrlpts,METHOD,V_PLANESPEED);


function CalculateManualValues(handles,h_ctrlpts,METHOD,V_PLANESPEED)
% This function updates the manual line and also calls the time calculation
% routine.

p0 = handles.p0;
p0(2:end-1,:) = cell2mat(get(nonzeros(h_ctrlpts),{'xdata','ydata'}));
XYdata = PointsToCurve(p0,METHOD,handles.IntegrationFineness);
set(handles.lineInteractiveCurve,'Xdata',XYdata(:,1),'Ydata',XYdata(:,2));
[t_interactive dist_interactive] = ...
calcTime(XYdata,handles.Ix,handles.Iy,V_PLANESPEED,METHOD,'interactive',handles);
updateManualText(handles,dist_interactive/t_interactive,dist_interactive,t_interactive);
drawnow;

% --- Executes on button press in pushbutton_Optimize.
function pushbutton_Optimize_Callback(hObject, eventdata, handles)
% This function carries out a global optimization using MultiStart.
%

% hObject    handle to pushbutton_Optimize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

METHOD = get(handles.popupmenu_METHOD,'value');
METHOD = getMETHODS(METHOD);
V_PLANESPEED = str2double(get(handles.edit_PLANESPEED,'string'));

% Set up the optimization problem:
opts = optimset('fmincon');
opts.Display = 'off';
opts.UseParallel = 'never';

Ix = handles.Ix;
Iy = handles.Iy;
SZX = handles.SZX;
SZY = handles.SZY;
npts = str2double(get(handles.edit_NUMPOINTS,'string'));

% These are nominal initial conditions (straight line path)
p0 = linspace(0,SZX,npts+2);
p0 = p0(2:end-1);
p0 = [p0; 0*p0+SZY/2];
p = [0 SZY/2; reshape(p0, 2,[])'; SZX SZY/2];

p0 = reshape(transpose(p(2:end-1,:)),1,[]); % Make p0 a vector

options = optimset('fmincon');
options.MaxFunEvals = npts*500;
options.Display = 'off';
options.UseParallel = 'never';
problem = createOptimProblem('fmincon','objective', ...
          @(X)calcTime(X,Ix,Iy,V_PLANESPEED,METHOD,'optimizer',handles), 'x0', p0(:),...
          'lb', 0.001+p0*0, ...
          'ub', repmat([SZX; SZY],size(p,1),1),...
          'options',options);
ms = MultiStart;
ms.Display = 'iter';
ms.UseParallel = 'always';
warning off;
try % Sometimes annoying warnings show up during the optimization...
    spmd
        warning off;
    end
end;
set(handles.pushbutton_Optimize,'String','Optimizing...','FontSize',8);
set(findall(gcbf,'type','uicontrol'),'enable','off'); % Disable UIcontrols
drawnow;
[popt, topt,EXITFLAG,OUTPUT,SOLUTIONS] = run(ms, problem,10); % Call MultiStart
set(findall(gcbf,'type','uicontrol'),'enable','on');
set(handles.pushbutton_Optimize,'String','Global');

p = [0 SZY/2; reshape(popt, 2,[])'; SZX SZY/2]; % reshape the solution into Nx2
hold all;
IntegrationFineness = handles.IntegrationFineness;
% Plot the optimal curve and update the "optimal" info box
XYdata = PointsToCurve(p,METHOD,IntegrationFineness);
lineOptimalCurve = plot(XYdata(:,1),XYdata(:,2),'k','linewidth',2);
handles.lineOptimalCurve = lineOptimalCurve;
set(lineOptimalCurve,'hittest','off');
for n = 2:(size(p,1)-1)
    hp(n) = plot(p(n,1),p(n,2),'ko');
end
set(hp(2:end),'tag','OptimalLine','hittest','off');

[t_opt dist_opt] = ...
calcTime(XYdata,handles.Ix,handles.Iy,V_PLANESPEED,METHOD,'interactive',handles);
updateOptText(handles,dist_opt/t_opt,dist_opt,t_opt);
%uistack([nonzeros(hp); lineOptimalCurve],'bottom')
%uistack([nonzeros(hp); lineOptimalCurve],'up')


function editSeed_Callback(hObject, eventdata, handles)
% hObject    handle to editSeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSeed as text
%        str2double(get(hObject,'String')) returns contents of editSeed as a double


% --- Executes during object creation, after setting all properties.
function editSeed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_PLANESPEED_Callback(hObject, eventdata, handles)
% When the Plane Airspeed is changed, this updates the plot.

% hObject    handle to edit_PLANESPEED (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PLANESPEED as text
%        str2double(get(hObject,'String')) returns contents of edit_PLANESPEED as a double
METHOD = get(handles.popupmenu_METHOD,'value');
METHOD = getMETHODS(METHOD);
V_PLANESPEED = str2double(get(handles.edit_PLANESPEED,'string'));
h_ctrlpts = handles.h_ctrlpts;

CalculateManualValues(handles,h_ctrlpts,METHOD,V_PLANESPEED)
p = handles.p0(2:end-1,:)';
t_straightLine = calcTime(p(:),handles.Ix,handles.Iy,V_PLANESPEED,METHOD,'optimizer',handles);
title(handles.AxesMain,sprintf('Straight Line Time: %d hours, %.1f minutes\n',floor(t_straightLine),rem(t_straightLine,1)*60));

% --- Executes during object creation, after setting all properties.
function edit_PLANESPEED_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PLANESPEED (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_METHOD.
function popupmenu_METHOD_Callback(hObject, eventdata, handles)
% When the Interpolation method is changed, this updates the plot.

% hObject    handle to popupmenu_METHOD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% ヒント: contents = cellstr(get(hObject,'String')) はセル配列として popupmenu_METHOD の内容を返します?B
%        contents{get(hObject,'Value')} は popupmenu_METHOD から選択した??ﾚを返します?B

METHOD = get(handles.popupmenu_METHOD,'value');
METHOD = getMETHODS(METHOD);
V_PLANESPEED = str2double(get(handles.edit_PLANESPEED,'string'));
h_ctrlpts = handles.h_ctrlpts;

CalculateManualValues(handles,h_ctrlpts,METHOD,V_PLANESPEED)


% --- Executes during object creation, after setting all properties.
function popupmenu_METHOD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_METHOD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_NUMPOINTS_Callback(hObject, eventdata, handles)
% When the number of points are changed, reset the interactive line

% hObject    handle to edit_NUMPOINTS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_NUMPOINTS as text
%        str2double(get(hObject,'String')) returns contents of edit_NUMPOINTS as a double
try
    makeInteractiveLine(handles);
catch %If there is bad input, reset it to 5 control points.
    set(hObject,'string','5');
    makeInteractiveLine(handles);
end

% --- Executes during object creation, after setting all properties.
function edit_NUMPOINTS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_NUMPOINTS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function varargout = calcTime(X,Ix,Iy,V_PLANESPEED,METHOD,caller,handles)
% This is the main function that actually calculates the line integral
% along the path.


p0 = handles.p0;
SZX = handles.SZX;
SZY = handles.SZY;
IntegrationFineness = handles.IntegrationFineness;

% If we are called from the optimization routine (caller = 'optimizer')
% then we need to interpolate the fine path from the input control points.
if strcmp(caller,'optimizer')
    p = [p0(1,:); reshape(X, 2,[])'; p0(end,:)];
    nP = size(p,1);
    plist = [interp1(1:nP,p(:,1),linspace(1,nP,IntegrationFineness)',METHOD,'extrap')...
        interp1(1:nP,p(:,2),linspace(1,nP,IntegrationFineness)',METHOD,'extrap')];
else %If we are in interactive mode, then we already have the fine point list
    plist = X;
end
%Do not let the path go out of the box
plist(:,1) = min(plist(:,1),SZX);
plist(:,1) = max(plist(:,1),0);
plist(:,2) = min(plist(:,2),SZY);
plist(:,2) = max(plist(:,2),0);

dP = diff(plist);

% We use a faster version of interp2 to interpolate the wind vector field
% at all the points in plist.
%Vwind_path = [interp2(Ix,plist(1:end-1,1)+1,plist(1:end-1,2)+1,'linear') interp2(Iy,plist(1:end-1,1)+1,plist(1:end-1,2)+1,'*linear')];
Vwind_path = [faster2dinterp(Ix,plist(1:end-1,:)+1) faster2dinterp(Iy,plist(1:end-1,:)+1)];

% Dot product the wind (Vwind_path) with the direction vector (dP) to get
% the tailwind/headwind contribution
V_add = (sum(Vwind_path.*dP,2))./sqrt(sum(dP.^2,2));

dx = sqrt(sum(dP.^2,2)); %dx is the length of each subinterval in plist
dT = dx./(V_PLANESPEED+V_add);  %dT = dP/dV
totalTime = 100*sum(dT);
%If the windspeed exceeds the airplane's airspeed, then we will get
%meaningless negative values.
if any(dT <= 0); 
    totalTime = Inf;
end
varargout{1} = totalTime;

% The optimizer only has one output, and will not need to output the
% path distance as we do here
if nargout > 1
    totalDist = sum(dx);
    varargout{2} = totalDist;
    drawnow;
end

function METHOD = getMETHODS(METHOD)
% List of interpolation methods for the popup menu.
list = {'linear', 'spline', 'pchip'};
METHOD = list{METHOD};

function updateManualText(handles,spd,d,t)
% Updates the text fields for the Manual path line.
if isempty(spd)
    set(handles.textManualSpeed,'String',sprintf('Average Speed: '))
    set(handles.textManualDistance,'String',sprintf('Distance: '))
    set(handles.textManualTime,'String',sprintf('Time: '));
else
    set(handles.textManualSpeed,'String',sprintf('Average Speed: %.3f kph',100*spd))
    set(handles.textManualDistance,'String',sprintf('Distance: %.3f km',100*d))
    set(handles.textManualTime,'String',sprintf('Time: %d hours, %.1f minutes\n',...
    floor(t),rem(t,1)*60)); %convert Hours --> Hours and Minutes...
end
drawnow;

function updateOptText(handles,spd,d,t)
% Updates the text fields for the Optimal path line.
if isempty(spd)
    set(handles.textOptSpeed,'String',sprintf('Average Speed: '))
    set(handles.textOptDistance,'String',sprintf('Distance: '))
    set(handles.textOptTime,'String',sprintf('Time: '));
else
    set(handles.textOptSpeed,'String',sprintf('Average Speed: %.3f kph',100*spd))
set(handles.textOptDistance,'String',sprintf('Distance: %.3f km',100*d))
set(handles.textOptTime,'String',sprintf('Time: %d hours, %.1f minutes\n',...
    floor(t),rem(t,1)*60));
end
drawnow;

function Iwind = makeRandomWindField(windFineness,SZX,SZY)
% This attempts to make a random vector field of wind

N = 50; % Various parameters used in generating a random "smooth" matrix
NL = 40; 
NP = 500;
rx = randn(NL,N);
rx = interpft(rx,NP);
ry = randn(NL,N);
ry = interpft(ry,NP);
I = (rx*ry');

[xgi,ygi] = meshgrid(linspace(1,2 + 498*windFineness,SZX+1),linspace(1,2 + 498*windFineness,SZY+1));
Iwind = interp2(1:500,1:500,I,xgi,ygi);



function edit_WindRoughness_Callback(hObject, eventdata, handles)
% hObject    handle to edit_WindRoughness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_WindRoughness as text
%        str2double(get(hObject,'String')) returns contents of edit_WindRoughness as a double


% --- Executes during object creation, after setting all properties.
function edit_WindRoughness_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_WindRoughness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_FieldX_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FieldX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_FieldX as text
%        str2double(get(hObject,'String')) returns contents of edit_FieldX as a double


% --- Executes during object creation, after setting all properties.
function edit_FieldX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FieldX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_FieldY_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FieldY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_FieldY as text
%        str2double(get(hObject,'String')) returns contents of edit_FieldY as a double


% --- Executes during object creation, after setting all properties.
function edit_FieldY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FieldY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function OUT = faster2dinterp(R,p)
% A helper function to do faster 2d interpolation on the wind vector field
p = p(:,[2 1]);
Aloc = floor(p);
localXY = p - Aloc;
Bloc = Aloc;
Bloc(:,1) = min(Bloc(:,1)+1,size(R,1));
Cloc = Aloc;
Cloc(:,2) = min(Cloc(:,2)+1,size(R,2));
Dloc = [min(Aloc(:,1),size(R,1)-1) min(Aloc(:,2),size(R,2)-1)]+1;

W1 = R(Aloc(:,1) + size(R,1)*Aloc(:,2)-size(R,1));
W2 = R(Bloc(:,1) + size(R,1)*Bloc(:,2)-size(R,1))-W1;
W3 = R(Cloc(:,1) + size(R,1)*Cloc(:,2)-size(R,1))-W1;
W4 = R(Dloc(:,1) + size(R,1)*Dloc(:,2)-size(R,1))-W1-W2-W3;
OUT = W1 + W2 .*localXY(:,1) + W3.*localXY(:,2) + W4.*localXY(:,1).*localXY(:,2);


% --- Executes on button press in pushbutton_OptimizeLocal.
function pushbutton_OptimizeLocal_Callback(hObject, eventdata, handles)
% This function carries out a local optimization with FMINCON starting at
% the current location of the control points for the interactive line.

% hObject    handle to pushbutton_OptimizeLocal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
METHOD = get(handles.popupmenu_METHOD,'value');
METHOD = getMETHODS(METHOD);
V_PLANESPEED = str2double(get(handles.edit_PLANESPEED,'string'));

opts = optimset('fmincon');
opts.Display = 'off';
opts.UseParallel = 'never';

Ix = handles.Ix;
Iy = handles.Iy;
SZX = handles.SZX;
SZY = handles.SZY;
npts = str2double(get(handles.edit_NUMPOINTS,'string'));

% Get the X and Y data from the control points.
p0 = cell2mat(get(nonzeros(handles.h_ctrlpts),{'xdata','ydata'}));

% end
options = optimset('fmincon');
options.MaxFunEvals = npts*5000;
options.Display = 'iter';
options.UseParallel = 'never';
options.OutputFcn = @(x, optimValues, state) outfun(x, optimValues, state,METHOD,handles,SZX,SZY);
options.TolFun = 5e-6;
options.TolX = 5e-6;

p0 = reshape(transpose(p0),1,[])'; % Make it a column vector

% Call FMINCON with calcTime as the objective function
[popt, topt] = fmincon(@(X)calcTime(X,Ix,Iy,V_PLANESPEED,METHOD,'optimizer',handles),p0(:),[],[],[],[],...
    p0*0,repmat([SZX; SZY],size(p0,1)/2,1),[],options);

p = [0 SZY/2; reshape(popt, 2,[])'; SZX SZY/2]; %Reshape back to Nx2
hold all;
% Draw the optimal line and update the "optimal" info box
IntegrationFineness = handles.IntegrationFineness;
XYdata = PointsToCurve(p,METHOD,IntegrationFineness);
lineOptimalCurve = plot(XYdata(:,1),XYdata(:,2),'k','linewidth',2);
handles.lineOptimalCurve = lineOptimalCurve;
set(lineOptimalCurve,'hittest','off');
for n = 2:(size(p,1)-1)
    hp(n) = plot(p(n,1),p(n,2),'ko');
end
set(hp(2:end-1),'tag','OptimalLine','hittest','off');

[t_opt dist_opt] = ...
calcTime(XYdata,handles.Ix,handles.Iy,V_PLANESPEED,METHOD,'interactive',handles);
updateOptText(handles,dist_opt/t_opt,dist_opt,t_opt);

function stop = outfun(x, optimValues, state,METHOD,handles,SZX,SZY)
% This is for interactive viewing of the optimization process. This will be
% called to update the graph when doing local optimization. This is not
% used during global optimization.

% Only run this function once every 10 times, so we don't slow down the
% optimization too much
persistent ctr
stop = false;
if isempty(ctr)
    ctr = 0;
end
ctr = ctr+1;
if mod(ctr,10)~= 1;
    return 
end
stop = false;

% Redraw the optimization solver's line based on the current point
popt = x;
p = [0 SZY/2; reshape(popt, 2,[])'; SZX SZY/2];
delete(findall(0,'tag','OptimalLine'));
hold all;
IntegrationFineness = handles.IntegrationFineness;
XYdata = PointsToCurve(p,METHOD,IntegrationFineness);
lineOptimalCurve = plot(XYdata(:,1),XYdata(:,2),'k','linewidth',2,'tag','OptimalLine');
handles.lineOptimalCurve = lineOptimalCurve;
set(lineOptimalCurve,'hittest','off');
set(0,'showhiddenhandles','on');
for n = 2:(size(p,1)-1)
    hp(n) = plot(p(n,1),p(n,2),'ko');
end
set(hp(2:end),'tag','OptimalLine','hittest','off');

XYdata(:,1) = min(XYdata(:,1),SZX);
XYdata(:,1) = max(XYdata(:,1),0);
XYdata(:,2) = min(XYdata(:,2),SZY);
XYdata(:,2) = max(XYdata(:,2),0);
dP = diff(XYdata);
dist_opt = sum(sqrt(sum(dP.^2,2)));
t_opt = optimValues.fval;
updateOptText(handles,dist_opt/t_opt,dist_opt,t_opt);
drawnow;
