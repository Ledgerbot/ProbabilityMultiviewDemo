%% DriveWallE
% This script demonstrates how you can move a visualization around a scene
% using MATLAB graphics objects and rigid body transformations.
%
%   M. Kutzer, 21Oct2019, USNA
clear all
close all
clc

%% Create video parameters
makeVideo = false;

if makeVideo
    vid(1) = VideoWriter('DriveWallE_3rdPerson_FixedOrient.mp4','MPEG-4');
    vid(2) = VideoWriter('DriveWallE_3rdPerson_2dOrient.mp4','MPEG-4');
    vid(3) = VideoWriter('DriveWallE_3rdPerson_3dUprightOrient.mp4','MPEG-4');
    vid(4) = VideoWriter('DriveWallE_3rdPerson_3dOrient.mp4','MPEG-4');
    for i = 1:numel(vid)
        open(vid(i));
    end
end

%% Create figure, axes, etc.
% Create figure
fig = figure('Name','Drive Wall-E');
% Set figure to 5.4" high x 12.67" (so we can insert in PPT w/o stretching)
set(fig,'Units','Inches','Position',[0.25,0.65,12.67,5.40]);
% Set figure background to white (so it will look good in PPT)
set(fig,'Color',[1,1,1]);

% Create visualization axes
axsVIS = axes('Parent',fig);
hold(axsVIS,'on');
daspect(axsVIS,[1 1 1]);
view(axsVIS,3);
xlim(axsVIS,[ -5.5,  5.5]);
ylim(axsVIS,[ -0.5,  5.5]);
zlim(axsVIS,[ -1.5,  3.0]);
set(axsVIS,'Position',[-0.1,-0.06,0.80,1.20]);
set(axsVIS,'Visible','off');

% Create text axes
axsMAT = axes('Parent',fig);
hold(axsMAT,'on');
daspect(axsMAT,[1 1 1]);
view(axsMAT,2);
xlim(axsMAT,[ -1.0,  1.0]);
ylim(axsMAT,[ -1.0,  1.0]);
set(axsMAT,'Position',[0.57,0.3,0.40,0.40]);
set(axsMAT,'Visible','off');
H = eye(4);
str = sprintf([...
    '%7.3f & %7.3f & %7.3f & %7.3f \\\\',...
    '%7.3f & %7.3f & %7.3f & %7.3f \\\\',...
    '%7.3f & %7.3f & %7.3f & %7.3f \\\\',...
    '0 & 0 & 0 & 1 \\\\'],reshape(H(1:3,:).',1,[]));

txt = text(0,0,['$H_{b}^{0} = \left( \begin{array}{rrrr} ',str,' \end{array} \right)$'],...
    'Parent',axsMAT,'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',20);

%set(axsMAT,'Visible','off');

% Create reference frames
frm_W = triad('Scale',1.00,'LineWidth',2,'AxisLabels',{'x_0','y_0','z_0'},'Parent',axsVIS);
frm_B = triad('Scale',0.75,'LineWidth',2,'AxisLabels',{'x_b','y_b','z_b'});
set(frm_B,'Parent',frm_W);

% Add light
lgt = addSingleLight(axsVIS);

%% Load Wall-E visualization
open('Wall-E.fig');
drawnow;
figTMP = findobj('Type','Figure','Name','Wall-E');
ptc = findobj('Type','Patch','Tag','Wall-E');
set(ptc,'Parent',frm_B);
delete(figTMP);

%% Create path
t = linspace(0,2*pi,300);
t = t - pi/2*ones(size(t));
% Position
x = 5.0*cos(1*t) + 0.0*ones(size(t));
y = 2.5*sin(1*t) + 2.5*ones(size(t));
z = 1.5*sin(4*t) + 0.5*ones(size(t));
% Velocity
dx = -5.0*sin(1*t);
dy =  2.5*cos(1*t);
dz =  6.0*cos(4*t);
% Acceleration
ddx =  -5.0*cos(1*t);
ddy =  -2.5*sin(1*t);
ddz = -24.0*sin(4*t);

% Plot path
plt = plot3(x,y,z,'m','LineWidth',1.5,'Parent',frm_W);

%% Have Wall-E follow the path

%% (1) with orientation fixed orientation
for i = 1:numel(t)
    H = eye(4);
    H(1:3,4) = [x(i); y(i); z(i)];
    set(frm_B,'Matrix',H);
    
    str = sprintf([...
        '%7.3f & %7.3f & %7.3f & %7.3f \\\\',...
        '%7.3f & %7.3f & %7.3f & %7.3f \\\\',...
        '%7.3f & %7.3f & %7.3f & %7.3f \\\\',...
        '0 & 0 & 0 & 1 \\\\'],reshape(H(1:3,:).',1,[]));
    set(txt,'String',['$H_{b}^{0} = \left( \begin{array}{rrrr} ',str,' \end{array} \right)$']);
    
    drawnow;
    if makeVideo
        frm = getframe(fig);
        writeVideo(vid(1),frm);
    end
    
end

%% (2) with orientation about a single axis
for i = 1:numel(t)
    H = eye(4);
    H(1:3,4) = [x(i); y(i); z(i)];
    z_hat = [0; 0; 1];
    x_hat = [dx(i); dy(i); 0];
    x_hat = x_hat./norm(x_hat);
    y_hat = cross(z_hat,x_hat);
    H(1:3,1:3) = [x_hat,y_hat,z_hat];
    set(frm_B,'Matrix',H);
    
    str = sprintf([...
        '%7.3f & %7.3f & %7.3f & %7.3f \\\\',...
        '%7.3f & %7.3f & %7.3f & %7.3f \\\\',...
        '%7.3f & %7.3f & %7.3f & %7.3f \\\\',...
        '0 & 0 & 0 & 1 \\\\'],reshape(H(1:3,:).',1,[]));
    set(txt,'String',['$H_{b}^{0} = \left( \begin{array}{rrrr} ',str,' \end{array} \right)$']);
    
    drawnow;
    if makeVideo
        frm = getframe(fig);
        writeVideo(vid(2),frm);
    end
    
end

%% (3) upright 3D orientation
for i = 1:numel(t)
    H = eye(4);
    H(1:3,4) = [x(i); y(i); z(i)];
    %z_hat = [ddx(i); ddy(i); ddz(i)];
    z_hat = [0; 0; 1];
    z_hat = z_hat./norm(z_hat);
    x_hat = [dx(i); dy(i); dz(i)];
    x_hat = x_hat./norm(x_hat);
    y_hat = cross(z_hat,x_hat);
    y_hat = y_hat./norm(y_hat);
    z_hat = cross(x_hat,y_hat);
    
    H(1:3,1:3) = [x_hat,y_hat,z_hat];
    set(frm_B,'Matrix',H);
    
    str = sprintf([...
        '%7.3f & %7.3f & %7.3f & %7.3f \\\\',...
        '%7.3f & %7.3f & %7.3f & %7.3f \\\\',...
        '%7.3f & %7.3f & %7.3f & %7.3f \\\\',...
        '0 & 0 & 0 & 1 \\\\'],reshape(H(1:3,:).',1,[]));
    set(txt,'String',['$H_{b}^{0} = \left( \begin{array}{rrrr} ',str,' \end{array} \right)$']);
    
    fprintf('%6.3f\n',det(H(1:3,1:3)));
    
    drawnow;
    if makeVideo
        frm = getframe(fig);
        writeVideo(vid(3),frm);
    end
    
end

%% (4) w/ 3D orientation
for i = 1:numel(t)
    H = eye(4);
    H(1:3,4) = [x(i); y(i); z(i)];
    z_hat = [ddx(i); ddy(i); ddz(i)];
    %z_hat = [0; 0; 1];
    z_hat = z_hat./norm(z_hat);
    x_hat = [dx(i); dy(i); dz(i)];
    x_hat = x_hat./norm(x_hat);
    y_hat = cross(z_hat,x_hat);
    y_hat = y_hat./norm(y_hat);
    z_hat = cross(x_hat,y_hat);
    
    H(1:3,1:3) = [x_hat,y_hat,z_hat];
    set(frm_B,'Matrix',H);
    
    str = sprintf([...
        '%7.3f & %7.3f & %7.3f & %7.3f \\\\',...
        '%7.3f & %7.3f & %7.3f & %7.3f \\\\',...
        '%7.3f & %7.3f & %7.3f & %7.3f \\\\',...
        '0 & 0 & 0 & 1 \\\\'],reshape(H(1:3,:).',1,[]));
    set(txt,'String',['$H_{b}^{0} = \left( \begin{array}{rrrr} ',str,' \end{array} \right)$']);
    
    fprintf('%6.3f\n',det(H(1:3,1:3)));
    
    drawnow;
    if makeVideo
        frm = getframe(fig);
        writeVideo(vid(4),frm);
    end
    
end
%% Close and delete videos
if makeVideo
    for i = 1:numel(vid)
        close(vid(i));
        delete(vid(i));
    end
end