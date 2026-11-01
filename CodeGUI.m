function varargout = CodeGUI(varargin)
%CODEGUI GUIDE-compatible front end for the shared frictional contact model.
% The GUI layout remains in Part 2/CodeGUI.fig; this controller is kept in
% the project root so it can run together with the main model files.
root = fileparts(mfilename('fullpath'));
guiState = struct('gui_Name',fullfile(root,'Part 2','CodeGUI'), ...
    'gui_Singleton',1,'gui_OpeningFcn',@CodeGUI_OpeningFcn, ...
    'gui_OutputFcn',@CodeGUI_OutputFcn,'gui_LayoutFcn',[],'gui_Callback',[]);
if nargin && ischar(varargin{1}), guiState.gui_Callback = str2func(varargin{1}); end
if nargout, [varargout{1:nargout}] = gui_mainfcn(guiState,varargin{:});
else, gui_mainfcn(guiState,varargin{:}); end
end

function CodeGUI_OpeningFcn(hObject,~,handles,~)
handles.output = hObject;
root = fileparts(mfilename('fullpath'));
addpath(root);
guidata(hObject,handles);
end

function varargout = CodeGUI_OutputFcn(~,~,handles), varargout{1} = handles.output; end

function start_button_Callback(~,~,handles)
r = simulate_from_gui(handles);
plot(r.time,r.modalDisplacement); grid on;
xlabel('Time (s)'); ylabel('Modal displacement (m)');
title(sprintf('Contact response, k_c = %.3g N/m',r.contactStiffness));
end

function pushbutton2_Callback(~,~,handles)
r = simulate_from_gui(handles);
plot(r.time,r.modalVelocity); grid on;
xlabel('Time (s)'); ylabel('Modal velocity (m/s)');
title(sprintf('Contact response, dissipated energy = %.4g J',r.energyDissipated));
end

function pushbutton3_Callback(~,~,handles)
p = gui_parameters(handles);
frequencies = linspace(0.2,2,80)*p.naturalFrequency;
amplitude = zeros(size(frequencies));
for i = 1:numel(frequencies)
    p.driveFrequency = frequencies(i);
    r = frictional_contact_response(p);
    nTime = numel(r.time);
    amplitude(i) = max(abs(r.modalDisplacement(round(nTime/2):nTime)));
end
plot(frequencies/(2*pi),amplitude); grid on;
xlabel('Drive frequency (Hz)'); ylabel('Steady amplitude (m)');
title('Frictional-contact frequency response');
end

function Graphics_Callback(~,~,handles)
%GRAPHICS_Callback Animate contact motion and update response plots together.
r = simulate_from_gui(handles);
beam = cantilever_modal_properties(210e9,7850,0.05*0.03, ...
    0.05*0.03^3/12,1);
frequency = linspace(0.2,2,50)*gui_parameters(handles).naturalFrequency;
frequencyAmplitude = frequency_response(gui_parameters(handles),frequency);

fig = figure('Name','Synchronized frictional-contact graphics', ...
    'NumberTitle','off','Color','w');
axAnimation = subplot(3,2,[1 2]);
axis(axAnimation,[-0.05 1.2 -1.2 1.2]); grid(axAnimation,'on');
xlabel(axAnimation,'Horizontal beam axis'); ylabel(axAnimation,'Scaled transverse displacement');
title(axAnimation,'Horizontal cantilever rubbing against a vertical surface'); hold(axAnimation,'on');
surfaceX = 1;
plot(axAnimation,[surfaceX surfaceX],[-1.2 1.2],'k-','LineWidth',2);
beamLine = plot(axAnimation,nan,nan,'b-','LineWidth',2);
contactBody = rectangle(axAnimation,'Position',[0.99 -0.12 0.12 0.12], ...
    'FaceColor','none','EdgeColor',[0.85 0.25 0.15],'LineWidth',1.5);
contactMarker = plot(axAnimation,nan,nan,'ko','MarkerFaceColor','none');

axDisp = subplot(3,2,3); plot(axDisp,r.time,r.modalDisplacement,'b'); grid(axDisp,'on');
xlabel(axDisp,'Time (s)'); ylabel(axDisp,'Displacement (m)'); title(axDisp,'Displacement'); hold(axDisp,'on');
dispMarker = plot(axDisp,nan,nan,'ro','MarkerFaceColor','none');
axHyst = subplot(3,2,4); plot(axHyst,r.contactDisplacement,r.contactForce,'b'); grid(axHyst,'on');
xlabel(axHyst,'Contact displacement (m)'); ylabel(axHyst,'Contact force (N)'); title(axHyst,'Contact hysteresis'); hold(axHyst,'on');
hystMarker = plot(axHyst,nan,nan,'ro','MarkerFaceColor','none');
axStress = subplot(3,2,5); plot(axStress,r.bendingStrain,r.bendingStress/1e6,'b'); grid(axStress,'on');
xlabel(axStress,'Bending strain'); ylabel(axStress,'Bending stress (MPa)'); title(axStress,'Stress-strain'); hold(axStress,'on');
stressMarker = plot(axStress,nan,nan,'ro','MarkerFaceColor','none');
axFreq = subplot(3,2,6); plot(axFreq,frequency/(2*pi),frequencyAmplitude,'b'); grid(axFreq,'on');
xlabel(axFreq,'Drive frequency (Hz)'); ylabel(axFreq,'Amplitude (m)'); title(axFreq,'Frequency response'); hold(axFreq,'on');

% Use the first-mode sine shape for the visual beam. The scale is visual only.
beamX = linspace(0,1,80);
beamShape = interp1(beam.x/beam.x(end),beam.modeShape,beamX);
visualScale = max([max(abs(r.contactDisplacement)), 1e-9]);
contactX = 1;
frameStep = max(1,ceil(numel(r.time)/600));
for k = 1:frameStep:numel(r.time)
    if ~ishandle(fig), break; end
    beamY = r.modalDisplacement(k)/visualScale*beamShape;
    contactY = r.contactDisplacement(k)/visualScale;
    set(beamLine,'XData',beamX,'YData',beamY);
    set(contactMarker,'XData',contactX,'YData',contactY);
    set(contactBody,'Position',[surfaceX-0.01,contactY-0.12,0.12,0.12]);
    set(dispMarker,'XData',r.time(k),'YData',r.modalDisplacement(k));
    set(hystMarker,'XData',r.contactDisplacement(k),'YData',r.contactForce(k));
    set(stressMarker,'XData',r.bendingStrain(k),'YData',r.bendingStress(k)/1e6);
    drawnow;
end
end

function amplitude = frequency_response(p,frequency)
amplitude = zeros(size(frequency));
for i = 1:numel(frequency)
    p.driveFrequency = frequency(i);
    response = frictional_contact_response(p);
    nTime = numel(response.time);
    amplitude(i) = max(abs(response.modalDisplacement(round(nTime/2):nTime)));
end
end

function pushbutton5_Callback(~,~,handles)
%PUSHBUTTON5_Callback Plot the bending stress-strain response.
r = simulate_from_gui(handles);
plot(r.bendingStrain,r.bendingStress/1e6,'LineWidth',1.1); grid on;
xlabel('Bending strain'); ylabel('Bending stress (MPa)');
title(sprintf('Bending stress-strain, E_{diss}=%.4g J',r.energyDissipated));
end

function result = simulate_from_gui(handles)
result = frictional_contact_response(gui_parameters(handles));
end

function p = gui_parameters(handles)
E = 210e9; rho = 7850; b = 0.05; h = 0.03; L = 1;
A = b*h; I = b*h^3/12;
beam = cantilever_modal_properties(E,rho,A,I,L);
p.mass = read_value(handles,'Mass',beam.modalMass);
p.damping = read_value(handles,'Damping',0.01);
p.duration = read_value(handles,'Time',2);
p.timeStep = read_value(handles,'Step',1e-3);
p.forceAmplitude = read_value(handles,'Amplitude',1);
p.driveFrequency = read_value(handles,'Omega',0.9*beam.naturalFrequency);
p.naturalFrequency = beam.naturalFrequency;
p.contactStiffness = 1e4;
p.gap = 0; p.frictionCoefficient = 0.35;
p.normalLoad = p.mass*9.81; p.modalShapeAtContact = 1;
p.frictionVelocity = 1e-5;
p.youngsModulus = E; p.beamLength = L; p.beamThickness = h;
p.curvatureAtStressLocation = beam.curvatureAtRootPerDisplacement;
end

function value = read_value(handles,name,defaultValue)
value = str2double(get(handles.(name),'String'));
if ~isfinite(value) || value <= 0, value = defaultValue; end
end

function Mass_Callback(~,~,~), end
function Time_Callback(~,~,~), end
function Step_Callback(~,~,~), end
function Damping_Callback(~,~,~), end
function Amplitude_Callback(~,~,~), end
function Omega_Callback(~,~,~), end
function Mass_CreateFcn(hObject,~,~), set_white(hObject); end
function Time_CreateFcn(hObject,~,~), set_white(hObject); end
function Step_CreateFcn(hObject,~,~), set_white(hObject); end
function Damping_CreateFcn(hObject,~,~), set_white(hObject); end
function Amplitude_CreateFcn(hObject,~,~), set_white(hObject); end
function Omega_CreateFcn(hObject,~,~), set_white(hObject); end
function set_white(hObject)
if ispc && isequal(get(hObject,'BackgroundColor'),get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
