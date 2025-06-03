%% ExampleCameraCalibration
% This script provides basic camera calibration tools
%
%   M. Kutzer, 03Mar2024, USNA
clear all
close all
clc

%% Initialize the camera
[cam,prv,handles] = initCamera;

camSettings = adjustCamera(cam);
%camSettings = adjustCamera(cam,camSettings,true);

%% Get Calibration Images
imBaseName = 'im';
calFolderName = 'ExampleCamCal';
nImages = 10;
getCalibrationImages(prv,imBaseName,calFolderName,nImages);

%% Camera Calibration
[cameraParams,imagesUsed] = ...
    calibrateCamera(imBaseName,calFolderName,nImages);

% NOTE: This can also be done using the MATLAB camera calibrator and
% manually defining the images used during calibration for reprojecting
%{
cameraCalibrator;
% Define Used Calibration Images
imagesUsed = [1,2,3,4,5,7,8,9,10];
%}

%% Parsing Camera Intrinsics and Extrinsics
% Parsing Intrinsic Matrix
A_c2m = cameraParams.IntrinsicMatrix.'; % <-- Note the transpose

% Parsing Extrinsic Matrices
n = cameraParams.NumPatterns; % <-- Total number of calibration images
for i = 1:n
    R_f2c = cameraParams.RotationMatrices(:,:,i).'; % <-- Note the transpose
    d_f2c = cameraParams.TranslationVectors(i,:).'; % <-- Note the transpose
    
    H_f2c{i} = [R_f2c, d_f2c; 0,0,0,1];	% <-- Each extrinsic matrix is contained in a cell
end

% Save parameters
save('ExampleCalibrationParameters.mat','A_c2m','H_f2c',...
    'calFolderName','cameraParams','imBaseName','imagesUsed','nImages');

%% Validating Your Calibration
% Recovering Body-fixed Fiducial Points
p_f = cameraParams.WorldPoints.'; % Parse x/y fiducial coordinates (note the transpose)
p_f(3,:) = 0;   % Fiducial z-coordinates (fiducial is 2D, so z is 0)
p_f(4,:) = 1;   % Make points homogeneous

% Reproject Body-fixed Fiducial Points
for i = 1:n
    % Define filename
    imName = sprintf('%s%03d.png',imBaseName,imagesUsed(i));
    % Load image
    im = imread( fullfile(calFolderName,imName) );
    % Plot image
    fig = figure('Name',imName);
    axs = axes('Parent',fig);
    img = imshow(im,'Parent',axs);
    hold(axs,'on');
    
    % Calculate projection matrix
    P_f2m = A_c2m * H_f2c{i}(1:3,:);
    % Project points using intrinsics and extrinsics
    tilde_p_m = P_f2m*p_f;	% <-- Scaled pixel coordinates 
    p_m = tilde_p_m./tilde_p_m(3,:); % <-- Pixel coordinates

    % Plot points
	% - Fiducial origin point
    plt0 = plot(axs,p_m(1,1),p_m(2,1),'ys','LineWidth',3,'MarkerSize',8); 
	% - All other fiducial points
    plti = plot(axs,p_m(1,2:end),p_m(2,2:end),'go','LineWidth',2,'MarkerSize',8);
    
    % Label points
    for j = 1:size(p_m,2)
        txt(j) = text(axs,p_m(1,j),p_m(2,j),sprintf('$p_{%d}^{m}$',j),...
            'Interpreter','Latex','Color','m','FontSize',14);
    end
end

%% Projecting 3D Computer Graphics into Images

% 1. Load the Wall-E visualization 
% Open Wall-E visualization and get figure handle
figWallE = open('Wall-E.fig');

% 2. Recover the Wall-E patch object 
% Recover the patch object from the Wall-E visualization figure handle
ptc_b = findobj(figWallE,'Type','patch','Tag','Wall-E');

% 3. Recover the vertices of Wall-E
p_b = ptc_b.Vertices.'; % <-- Note the transpose
p_b(4,:) = 1; % Make points homogeneous positions

% 4. Define Frame $o$ relative to the fiducial frame $f$
% Define checkerboard info
[boardSize,squareSize] = checkerboardPoints2boardSize( cameraParams.WorldPoints );
% Define Frame o relative to the fiducial frame f
x_o2f = squareSize*(boardSize(2)-2)/2; % x-offset of frame o wrt frame f
y_o2f = squareSize*(boardSize(1)-2)/2; % y-offset of frame o wrt frame f
H_o2f = Tx(x_o2f)*Ty(y_o2f)*Rx(pi);

% 5. Define Wall-E’s body-fixed frame relative to frame o. 
%    Note that we will initially use the identity, and update this frame  
%    when we make Wall-E drive.
H_b2o = eye(4);

% 6. Choose an image index to project Wall-E onto 
%    (e.g. i = 5 for the 5th image)
i = 5;

% 7. Define the projection matrix associated with this image:
% Define applicable extrinsics
H_b2c = H_f2c{i}*H_o2f*H_b2o;
% Projection matrix
P_b2m = A_c2m * H_b2c(1:3,:);

% 8. Project Wall-E's vertices onto the image
% Project points using intrinsics and extrinsics
tilde_p_m = P_b2m*p_b; % <-- Scaled pixel coordinates
p_m = tilde_p_m./tilde_p_m(3,:); % <-- Pixel coordinates

% 9. Load and plot image
% Define filename
imName = sprintf('%s%03d.png',imBaseName,imagesUsed(i));
% Load image
im = imread( fullfile(calFolderName,imName) );
% Plot image
fig = figure('Name',imName);
axs = axes('Parent',fig);
img = imshow(im,'Parent',axs);
hold(axs,'on');

% [DEBUG] Show fiducial and "body" frames ---------------------------------
sc = squareSize*1.5;
plt_b2m = projectTriad(axs,P_b2m,sc);
plt_f2m = projectTriad(axs,A_c2m * H_f2c{i}(1:3,:),sc);
% -------------------------------------------------------------------------

% 10. Render Wall-E's projection in the image
% Place Wall-E in the image
ptc_m = copyobj(ptc_b,axs);
ptc_m.Vertices = p_m(1:2,:).';

%% Improving the Visualization
% [ALLOWS RUNNING THIS SECTION ONLY]
% Open Wall-E visualization and get figure handle
if ~exist('figWallE','var') || ~ishandle(figWallE)
    figWallE = open('Wall-E.fig');
end

% 1. Recover the Wall-E patch object 
% Recover the patch object from the Wall-E visualization figure handle
ptc_b = findobj(figWallE,'Type','patch','Tag','Wall-E');

% 2. Recover the vertices of Wall-E
p_s = ptc_b.Vertices.'; % <-- Note the transpose
p_s(4,:) = 1; % Make points homogeneous positions

% 3. Define Frame o relative to the fiducial frame f
% Define checkerboard info
[boardSize,squareSize] = checkerboardPoints2boardSize( cameraParams.WorldPoints );
% Define Frame o relative to the fiducial frame f
x_o2f = squareSize*(boardSize(2)-2)/2; % x-offset of frame o wrt frame f
y_o2f = squareSize*(boardSize(1)-2)/2; % y-offset of frame o wrt frame f
H_o2f = Tx(x_o2f)*Ty(y_o2f)*Rx(pi);
% Scale Wall-E for square size
scWallE = 2;
p_b = Sx(scWallE*squareSize)*Sy(scWallE*squareSize)*Sz(scWallE*squareSize)*p_s;

% 4. Define Wall-E’s body-fixed frame relative to frame o. Note that we will initially use the identity, and
% update this frame when we make Wall-E drive.
H_b2o = eye(4);

% 5. Choose an image index to project Wall-E onto (e.g. i = 5 for the 5th image)
i = 5;

% 6. Define the projection matrix associated with this image:
% Define applicable extrinsics
H_b2c = H_f2c{i}*H_o2f*H_b2o;
% Projection matrix
P_b2m = A_c2m * H_b2c(1:3,:);

% 7. Project Wall-E’s vertices onto the image with false depth
% Project points with added false depth
p_m_falseDepth = projectWithFalseDepth(p_b,P_b2m);

% 8. Load and plot image
% Define filename
imName = sprintf('%s%03d.png',imBaseName,imagesUsed(i));
% Load image
im = imread( fullfile(calFolderName,imName) );
% Plot image
fig = figure('Name',imName);
axs = axes('Parent',fig);
img = imshow(im,'Parent',axs);
hold(axs,'on');

% 9. Add a light to the scene to better show Wall-E
addSingleLight(axs);

% [DEBUG] Show fiducial and "body" frames ---------------------------------
sc = squareSize*1.5;
plt_b2m = projectTriad(axs,P_b2m,sc);
plt_f2m = projectTriad(axs,A_c2m * H_f2c{i}(1:3,:),sc);
% -------------------------------------------------------------------------

% 10. Render Wall-E’s projection in the image
% Place Wall-E in the image
ptc_m = copyobj(ptc_b,axs);
ptc_m.Vertices = p_m_falseDepth.';

%% Animating the Improved Visualization
% [ALLOWS RUNNING THIS SECTION ONLY]
% Open Wall-E visualization and get figure handle
if ~exist('figWallE','var') || ~ishandle(figWallE)
    figWallE = open('Wall-E.fig');
end

% Define position information
d_b2o = []; x_hat = [];
r = squareSize*min(boardSize)/2;
phi = linspace(0,2*pi,50);
d_b2o(1,:) = r*cos(phi);  % x-position
d_b2o(2,:) = r*sin(phi);  % y-position
d_b2o(3,:) = 0;           % z-position 

% Define orientation information
z_hat = [0;0;1];        % z-direction
% x-direction
x_hat(1,:) = -sin(phi);
x_hat(2,:) =  cos(phi);
x_hat(3,:) =  0;

% Define rigid body transforms
H_b2o = cell(1,numel(phi));
for k = 1:numel(phi)
    % Define Wall-E's pose relative to the center frame
    y_hat = cross(z_hat,x_hat(:,k));
    % Define orientation
    R_b2o = [x_hat(:,k),y_hat,z_hat];
    % Define pose
    H_b2o{k} = [R_b2o, d_b2o(:,k); 0,0,0,1];
end

% Recover the patch object from the Wall-E visualization figure handle
ptc_b = findobj(figWallE,'Type','patch','Tag','Wall-E');

% Recover the vertices of Wall-E
p_s = ptc_b.Vertices.'; % <-- Note the transpose
p_s(4,:) = 1; % Make points homogeneous positions

% Define Frame o relative to the fiducial frame f
% Define checkerboard info
[boardSize,squareSize] = checkerboardPoints2boardSize( cameraParams.WorldPoints );
% Define Frame o relative to the fiducial frame f
x_o2f = squareSize*(boardSize(2)-2)/2; % x-offset of frame o wrt frame f
y_o2f = squareSize*(boardSize(1)-2)/2; % y-offset of frame o wrt frame f
H_o2f = Tx(x_o2f)*Ty(y_o2f)*Rx(pi);
% Scale Wall-E for square size
scWallE = 2;
p_b = Sx(scWallE*squareSize)*Sy(scWallE*squareSize)*Sz(scWallE*squareSize)*p_s;

% Recover the patch object from the Wall-E visualization figure handle
p_b(end,:) = []; % Do not make the points homogeous, this will make 
                 % projectWithFalseDepth run a little faster.

% Define image
i = 5;
% Define filename
imName = sprintf('%s%03d.png',imBaseName,imagesUsed(i));
% Load image
im = imread( fullfile(calFolderName,imName) );
% Plot image
fig = figure('Name',imName);
axs = axes('Parent',fig);
img = imshow(im,'Parent',axs);
hold(axs,'on');
addSingleLight(axs);
ptc_m = copyobj(ptc_b,axs);

% [DEBUG] Show fiducial and "body" frames ---------------------------------
sc = squareSize*1.5;
H_o2c = H_f2c{i}*H_o2f;
P_o2m = A_c2m*H_o2c(1:3,:);
plt_o2m = projectTriad(axs,P_o2m,sc);
plt_f2m = projectTriad(axs,A_c2m * H_f2c{i}(1:3,:),sc);
sc = squareSize*0.5;
for k = 1:numel(H_b2o)
    H_b2c = H_f2c{i}*H_o2f*H_b2o{k};
    P_b2m_k = A_c2m * H_b2c(1:3,:);
    projectTriad(axs,P_b2m_k,sc);
end
% -------------------------------------------------------------------------

for k = 1:numel(H_b2o)
	% Define applicable extrinsics
	H_b2c = H_f2c{i}*H_o2f*H_b2o{k};
	
	% Define projection
    P_b2m = A_c2m * H_b2c(1:3,:);
    
    % Project points with added false depth
    p_m_falseDepth = projectWithFalseDepth(p_b,P_b2m);
    
    % Update Wall-E
    ptc_m.Vertices = p_m_falseDepth.';
    drawnow
end

%% Animating the Improved Visualization (Creating a video)
% [ALLOWS RUNNING THIS SECTION ONLY]
% Open Wall-E visualization and get figure handle
if ~exist('figWallE','var') || ~ishandle(figWallE)
    figWallE = open('Wall-E.fig');
end

% Define position information
d_b2o = []; x_hat = [];
r = squareSize*min(boardSize)/2;
phi = linspace(0,2*pi,50);
d_b2o(1,:) = r*cos(phi);  % x-position
d_b2o(2,:) = r*sin(phi);  % y-position
d_b2o(3,:) = 0;           % z-position 

% Define orientation information
z_hat = [0;0;1];        % z-direction
% x-direction
x_hat(1,:) = -sin(phi);
x_hat(2,:) =  cos(phi);
x_hat(3,:) =  0;

% Define rigid body transforms
H_b2o = cell(1,numel(phi));
for k = 1:numel(phi)
    % Define Wall-E's pose relative to the center frame
    y_hat = cross(z_hat,x_hat(:,k));
    % Define orientation
    R_b2o = [x_hat(:,k),y_hat,z_hat];
    % Define pose
    H_b2o{k} = [R_b2o, d_b2o(:,k); 0,0,0,1];
end

% Recover the patch object from the Wall-E visualization figure handle
ptc_b = findobj(figWallE,'Type','patch','Tag','Wall-E');

% Define Frame o relative to the fiducial frame f
% Define checkerboard info
[boardSize,squareSize] = checkerboardPoints2boardSize( cameraParams.WorldPoints );
% Define Frame o relative to the fiducial frame f
x_o2f = squareSize*(boardSize(2)-2)/2; % x-offset of frame o wrt frame f
y_o2f = squareSize*(boardSize(1)-2)/2; % y-offset of frame o wrt frame f
H_o2f = Tx(x_o2f)*Ty(y_o2f)*Rx(pi);
% Scale Wall-E for square size
scWallE = 2;
p_b = Sx(scWallE*squareSize)*Sy(scWallE*squareSize)*Sz(scWallE*squareSize)*p_s;

% Recover the patch object from the Wall-E visualization figure handle
p_b(end,:) = []; % Do not make the points homogeous, this will make 
                 % projectWithFalseDepth run a little faster.

% Define image
i = 5;
% Define filename
imName = sprintf('%s%03d.png',imBaseName,imagesUsed(i));
% Load image
im = imread( fullfile(calFolderName,imName) );
% Plot image
fig = figure('Name',imName,'Color',[1 1 1]);
axs = axes('Parent',fig);
img = imshow(im,'Parent',axs);
hold(axs,'on');
addSingleLight(axs);
ptc_m = copyobj(ptc_b,axs);

% /////////////////////////////////////////////////////////////////////////
% NEW CODE
% Initialize an MPEG-4 video writer object
vid = VideoWriter('WallE_Driving.mp4','MPEG-4');
% Open video writer object
open(vid);
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

for k = 1:numel(H_b2o)
    % Define applicable extrinsics
    H_b2c = H_f2c{i}*H_o2f*H_b2o{k};
    
    % Define projection
    P_b2m = A_c2m * H_b2c(1:3,:);
    
    % Project points with added false depth
    p_m_falseDepth = projectWithFalseDepth(p_b,P_b2m);
    
    % Update Wall-E
    ptc_m.Vertices = p_m_falseDepth.';
    drawnow
    
    
    % /////////////////////////////////////////////////////////////////////
    % NEW CODE (figure handle)
    % "Take a picture" of the figure handle
    % (this assumes "fig" is the figure object you want in the video)
    frame = getframe(fig);
    % Write the frame to the video
    writeVideo(vid,frame);
    % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
    
    %{
    % /////////////////////////////////////////////////////////////////////
    % "Take a picture" of the axes handle
    % (this assumes "axs" is the axes object you want in the video)
    frame = getframe(axs);
    % Write the frame to the video
    writeVideo(vid,frame);
    % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    %}
end

% /////////////////////////////////////////////////////////////////////////
% NEW CODE
% Close video writer object
close(vid);
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

