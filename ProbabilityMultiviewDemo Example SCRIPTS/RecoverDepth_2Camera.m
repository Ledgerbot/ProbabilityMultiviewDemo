%% RecoverDepth_2Camera
% Simulated example of recovering depth of a segmented ball using two
% cameras
%
%   M. Kutzer, 23Mar2022, USNA

clear all
close all
clc

makeVideo = false;

if makeVideo
    % Setup to create video
    vidTitle = sprintf('RecoverDepth_2Camera_3D.mp4');
    vidObj_3D = VideoWriter(vidTitle,'MPEG-4');
    open(vidObj_3D);

    % Setup to create video
    vidTitle = sprintf('RecoverDepth_2Camera_cA.mp4');
    vidObj_c(1) = VideoWriter(vidTitle,'MPEG-4');
    open(vidObj_c(1));

    % Setup to create video
    vidTitle = sprintf('RecoverDepth_2Camera_cB.mp4');
    vidObj_c(2) = VideoWriter(vidTitle,'MPEG-4');
    open(vidObj_c(2));
end

%% Create figure, axes, etc.
fig = figure('Name','3D Simulation');
axs = axes('Parent',fig);
hold(axs,'on');
daspect(axs,[1 1 1]);

% Add lighting
lgt(1) = light(axs,'Style','Local','Position',[-100,-100,2000]);
lgt(2) = light(axs,'Style','Local','Position',[ 100,-100,2000]);
lgt(3) = light(axs,'Style','Local','Position',[ 100, 100,2000]);
lgt(4) = light(axs,'Style','Local','Position',[-100, 100,2000]);
set(lgt,'Color',[0.25,0.25,0.25]);

%% Adjust figure and axes
% Adjust view
view(axs,[0,25]);

% Adjust axes limits
xlim(axs,[-900, 900]);
ylim(axs,[-900, 900])
zlim(axs,[   0,1200]);

% Label axes
xlabel(axs,'x (mm)');
ylabel(axs,'y (mm)');
zlabel(axs,'z (mm)');

% Set figure position and color for video
set(fig,'Color',[1 1 1],...
    'Units','Inches','Position',[0.25,0.65,6.27,4.70]);

%% Load camera parameters
% Check version of MATLAB
v = matlabRelease;
d_v = datetime(v.Date);
d_s = datetime('31-Dec-2023');
if d_v < d_s
    % MATLAB R2023b or older
    load('ExampleCameraCalibration2023.mat');
else
    % MATLAB R2024a or newer
    load('ExampleCameraCalibration2024.mat');
end

% -> Pretend we have two sets of camera parameters
%    (we will calibrate both camera seperately to get this!)
tmp = cameraParams;
cameraParams = {};
cameraParams{1} = tmp;
cameraParams{2} = tmp;

% Create limited copy of cameraParams.mat
% -> Doing this ignores distortion when simulating images
% -> We are doing this because of a current limitation in simulateImage.m
%    that causes projected faces of patch objects to appear jagged along
%    the edge of images.
% -> We will use params for all simulated images (simulateImage.m)
for i = 1:numel(cameraParams)
    params{i}.IntrinsicMatrix = cameraParams{i}.IntrinsicMatrix;
    params{i}.ImageSize       = cameraParams{i}.ImageSize;
end

%% Place & visualize simulated cameras
H_c2a{1} = Tz(1100)*Tx( 200)*Ry(pi)*Rz(pi/2)*Rx( deg2rad(10));
h_c2a(1) = drawDFKCam(axs,H_c2a{1});

H_c2a{2} = Tz(1100)*Tx(-200)*Ry(pi)*Rz(pi/2)*Rx(-deg2rad(10));
h_c2a(2) = drawDFKCam(axs,H_c2a{2});

%% Visualize camera FOV
depth = 1300;
fov_c(1) = plotCameraFOV(h_c2a(1),cameraParams{1},depth);
set(fov_c(1),'FaceAlpha',0.1,'FaceColor','r');

fov_c(2) = plotCameraFOV(h_c2a(2),cameraParams{2},depth);
set(fov_c(2),'FaceAlpha',0.1,'FaceColor','b');

%% Simulate fiducial (we will a red ball)
% Define ball radius
% -> We do not need to know the diameter of the ball for this multiview
%    geometry approach!
r = 25; % mm

% Simulate ball
[p_f,h_f2a] = plotSphere(axs,[0,0,0,r]);

% Make ball red and opaque
set([p_f,p_f],'FaceColor','r','FaceAlpha',1);

% Hide ball frames
hideTriad(h_f2a);

% Move ball
H_f2a = Tz(100);
set(h_f2a,'Matrix',H_f2a);

%% Adjust axes limits
axis(axs,'tight');
zz = zlim(axs);
zlim(axs,[0,zz(2)]);

%% Take simulated images
H_a2c{1} = invSE(H_c2a{1});
H_a2c{2} = invSE(H_c2a{2});

% Hide the camera and FOV simulations
set(h_c2a,'Visible','off');
% Take simulated images
im{1} = simulateImage(axs,params{1},H_a2c{1});
im{2} = simulateImage(axs,params{2},H_a2c{2});
% Show the camera and FOV simulations
set(h_c2a,'Visible','on');

%% Show simulated images
for i = 1:2
    fig2D(i) = figure('Name',sprintf('Simulated Image %d',i));
    axs2D(i) = axes('Parent',fig2D(i));
    img(i) = imshow(im{i},'Parent',axs2D(i));
    ttl(i) = title(axs2D(i),sprintf('Camera %d',i));
end
hold(axs2D,'on');
set(axs2D,'Visible','on');
xlabel(axs2D,'x (pixels)');
ylabel(axs2D,'y (pixels)');

% Set figure position and color for video
set(fig2D(1),'Color',[1 1 1],...
    'Units','Inches','Position',[6.80,5.50,6.27,4.70]);
set(fig2D(2),'Color',[1 1 1],...
    'Units','Inches','Position',[0.25,5.50,6.27,4.70]);

%% Change ball alpha for visualization
set(p_f,'FaceAlpha',0.3);

%% Get ball position in images & visualize
for i = 1:2
    % Get ball position
    % -> Note that we are *not* using the intrinsic matrix or ball radius
    [p_m{i},~,ps_m(i)] = propsRedBall(im{i});

    % Visualize ball segmentation & centroid
    psp2D(i) = plot(ps_m(i),'Parent',axs2D(i),'FaceColor','m','FaceAlpha',0.5);
    plt2D(i) = plot(p_m{i}(1),p_m{i}(2),'+k','Parent',axs2D(i));
end

%% Define parameters to recover depth
% Parse intrinsic matrices
for i = 1:2
    A_c2m{i} = cameraParams{i}.IntrinsicMatrix.'; % <--- Note the transpose
end

% Define extrinsics
% -> This will be recovered experimentally using images of checkerboard
%    fiducials in a common pose
H_c22c1 = invSE(H_c2a{1})*H_c2a{2};
% -> Isolate rotation & translation
R_c22c1 = H_c22c1(1:3,1:3);
d_c22c1 = H_c22c1(1:3,4);

%% Recover depth
for i = 1:2
    tilde_p_c{i} = (A_c2m{i}^(-1))*p_m{i};
end
% Two-camera form
%z_c_i = pinv( [tilde_p_c{1}, -R_c22c1*tilde_p_c{2}] )*d_c22c1;
% General N-camera form
H_c12c2 = invSE(H_c22c1);
R_c12c2 = H_c12c2(1:3,1:3);
d_c12c2 = H_c12c2(1:3,4);
M = [...
    tilde_p_c{1}, -R_c22c1*tilde_p_c{2};...
    -R_c12c2*tilde_p_c{1}, tilde_p_c{2}];
z_c_i = pinv( M )*[d_c22c1; d_c12c2];

%% Visualize lines associated with depth estimate
colors = 'rb';
for i = 1:2
    % Interpolate "parameterization"
    z_c_lin(i,:) = linspace(0,depth,100);
    % Plot line
    % -> NOTE: tilde_p_c represents the coefficients for a line containing
    %          the origin and the point in 3D space (assuming zero error)!
    plt_lin(i) = plot3(h_c2a(i),...
        tilde_p_c{i}(1)*z_c_lin(i,:),...
        tilde_p_c{i}(2)*z_c_lin(i,:),...
        tilde_p_c{i}(3)*z_c_lin(i,:),...
        colors(i),'LineWidth',1);
end

%% Recover 3D position relative to each camera
for i = 1:2
    p_c_i{i} = z_c_i(i)*tilde_p_c{i};
end

%% Plot recovered 3D position relative to each camera
% Plot recovered points
plt_rec(1) = plot3(h_c2a(1),...
    p_c_i{1}(1),p_c_i{1}(2),p_c_i{1}(3),['x',colors(1)]);
plt_rec(2) = plot3(h_c2a(2),...
    p_c_i{2}(1),p_c_i{2}(2),p_c_i{2}(3),['+',colors(2)]);
% Plot segment connecting them
p_a_i(:,1) = H_c2a{1}*[p_c_i{1}; 1];
p_a_i(:,2) = H_c2a{2}*[p_c_i{2}; 1];
plt_seg = plot3(axs,...
    p_a_i(1,:),p_a_i(2,:),p_a_i(3,:),'k','LineWidth',1);

%% Define common p_c relative to each camera
p_c{1} = mean( [p_c_i{1}, H_c22c1(1:3,:)*[p_c_i{2}; 1]], 2 );
H_c12c2 = invSE(H_c22c1);
p_c{2} = H_c12c2(1:3,:)*[p_c{1}; 1];

%% Plot common 3D position relative to each camera
plt_com(1) = plot3(h_c2a(1),...
    p_c{1}(1),p_c{1}(2),p_c{1}(3),'ok');
plt_com(2) = plot3(h_c2a(2),...
    p_c{2}(1),p_c{2}(2),p_c{2}(3),'ok');

%% Project common 3D position into each image & plot
for i = 1:2
    tilde_p_m_com{i} = A_c2m{i}*p_c{i};
    p_m_com{i} = tilde_p_m_com{i}./tilde_p_m_com{i}(3,:);
    plt2D_com(i) = plot(axs2D(i),p_m_com{i}(1),p_m_com{i}(2),'ok');
end

%% Get actual ball position relative to each camera
for i = 1:2
    H_f2c{i} = invSE(H_c2a{i})*H_f2a;
    f_c{i} = H_f2c{i}(1:3,:)*[0;0;0;1];
end

%% Plot true 3D position relative to each camera
plt_tru(1) = plot3(h_c2a(1),...
    f_c{1}(1),f_c{1}(2),f_c{1}(3),'xk');
plt_tru(2) = plot3(h_c2a(2),...
    f_c{2}(1),f_c{2}(2),f_c{2}(3),'xk');

%% Project true 3D position into each image & plot
for i = 1:2
    tilde_p_m_tru{i} = A_c2m{i}*f_c{i};
    p_m_tru{i} = tilde_p_m_tru{i}./tilde_p_m_tru{i}(3,:);
    plt2D_tru(i) = plot(axs2D(i),p_m_tru{i}(1),p_m_tru{i}(2),'xk');
end

%% Create legends
lgnd = legend(axs,[plt_tru, plt_rec,plt_seg,plt_com],...
    'True Center (CamA ref)','True Center (CamB ref)',...
    'Recovered Center (CamA)','Recovered Center (CamB)',...
    'Segment Connecting Centers',...
    'Common Center (CamA ref)','Common Center (CamB ref)');

for i = 1:2
    lgnd2D(i) = legend(axs2D(i),[plt2D_tru(i), plt2D(i),plt2D_com(i)],...
        'True Center',...
        'Recovered Center',...
        'Common Center');
end

%% Zoom in on intersections
% 3D axes
% Define desired zoom (3D)
p_a(:,1) = H_c2a{1}*[p_c{1}; 1];
p_a(:,2) = H_c2a{2}*[p_c{2}; 1];
f_a(:,1) = H_c2a{1}*[f_c{1}; 1];
f_a(:,2) = H_c2a{2}*[f_c{2}; 1];
for i = 1:3
    limZ(i,:) = [...
        min([p_a_i(i,:),p_a(i,:),f_a(i,:)]),...
        max([p_a_i(i,:),p_a(i,:),f_a(i,:)])] + [-10,10];
end
% Define current zoom (3D)
limN(1,:) = xlim(axs);
limN(2,:) = ylim(axs);
limN(3,:) = zlim(axs);
% Interpolate zoom
n = 15*30; % <--- Zoom view frames
for i = 1:3
    limA{i}(1,:) = linspace(limN(i,1),limZ(i,1),n);
    limA{i}(2,:) = linspace(limN(i,2),limZ(i,2),n);
end

% 2D axes
for k = 1:2
    % Define desired zoom (2D)
    for i = 1:2
        lim2dZ{k}(i,:) = [...
            min([p_m_com{k}(i,:),p_m{k}(i,:),p_m_tru{k}(i,:)]),...
            max([p_m_com{k}(i,:),p_m{k}(i,:),p_m_tru{k}(i,:)])] + [-5,5];
    end
    % Define current zoom (2D)
    lim2dN{k}(1,:) = xlim(axs2D(k));
    lim2dN{k}(2,:) = ylim(axs2D(k));
    % Interpolate zoom
    for i = 1:2
        lim2dA{k}{i}(1,:) = linspace(lim2dN{k}(i,1),lim2dZ{k}(i,1),n);
        lim2dA{k}{i}(2,:) = linspace(lim2dN{k}(i,2),lim2dZ{k}(i,2),n);
    end
end

%% Take still video
for i = 1:30
    drawnow
    if makeVideo
        % Grab frame for video
        frame = getframe(fig);
        writeVideo(vidObj_3D,frame);
        for k = 1:2
            frame = getframe(fig2D(k));
            writeVideo(vidObj_c(k),frame);
        end
    end
end

% Zoom
for i = 1:n
    xlim(axs,[limA{1}(:,i).']);
    ylim(axs,[limA{2}(:,i).']);
    zlim(axs,[limA{3}(:,i).']);
    for k = 1:2
        xlim(axs2D(k),[lim2dA{k}{1}(:,i).']);
        ylim(axs2D(k),[lim2dA{k}{2}(:,i).']);
    end
    drawnow;
    if makeVideo
        % Grab frame for video
        frame = getframe(fig);
        writeVideo(vidObj_3D,frame);
        for k = 1:2
            frame = getframe(fig2D(k));
            writeVideo(vidObj_c(k),frame);
        end
    end
end

%% Rotate view
[az,el] = view(axs);
m = 15*30; % <--- Rotate view frames
azALL = linspace(az,az+360,m);
for i = 1:m
    view(axs,[azALL(i),el]);
    drawnow
    if makeVideo
        % Grab frame for video
        frame = getframe(fig);
        writeVideo(vidObj_3D,frame);
        for k = 1:2
            frame = getframe(fig2D(k));
            writeVideo(vidObj_c(k),frame);
        end
    end
end

if makeVideo
    % Close video obj
    close(vidObj_3D);
    close(vidObj_c(1));
    close(vidObj_c(2));
end