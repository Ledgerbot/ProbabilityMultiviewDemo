%% SCRIPT_Lab07_Solution
clear all
close all
clc

%% Initialize the camera
[cam,prv,handles] = initCamera;

camSettings = adjustCamera(cam);
%camSettings = adjustCamera(cam,camSettings,true);

%% Acquiring Calibration Images
imBaseName = 'im';
calFolderName = 'ExampleCamCal';
images = 10;
getCalibrationImages(prv,imBaseName,calFolderName,10);

%% Camera Calibration
cameraCalibrator;

%% Defining Used Calibration Images
imagesUsed = [1,2,3,4,5,7,8,9,10];

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
save('Lab07_calibrationParameters.mat','A_c2m','H_f2c',...
    'calFolderName','cameraParams','imBaseName','imagesUsed');

%% Validating Your Calibration
% [ALLOWS RUNNING THIS SECTION ONLY] vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
load('Lab07_calibrationParameters.mat');
n = cameraParams.NumPatterns; % <-- Total number of calibration images
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
% [ALLOWS RUNNING THIS SECTION ONLY] vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
load('Lab07_calibrationParameters.mat');
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

% 1. Load the Wall-E visualization 
% Open Wall-E visualization and get figure handle
figWallE = open('Lab07_Wall-E.fig');

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
% [ALLOWS RUNNING THIS SECTION ONLY] vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
load('Lab07_calibrationParameters.mat');
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

% [ALLOWS RUNNING THIS SECTION ONLY]
% Open Wall-E visualization and get figure handle
if ~exist('figWallE','var') || ~ishandle(figWallE)
    figWallE = open('Lab07_Wall-E.fig');
end

% 1. Recover the Wall-E patch object 
% Recover the patch object from the Wall-E visualization figure handle
ptc_b = findobj(figWallE,'Type','patch','Tag','Wall-E');

% 2. Recover the vertices of Wall-E
p_b = ptc_b.Vertices.'; % <-- Note the transpose
p_b(4,:) = 1; % Make points homogeneous positions

% 3. Define Frame o relative to the fiducial frame f
% Define checkerboard info
[boardSize,squareSize] = checkerboardPoints2boardSize( cameraParams.WorldPoints );
% Define Frame o relative to the fiducial frame f
x_o2f = squareSize*(boardSize(2)-2)/2; % x-offset of frame o wrt frame f
y_o2f = squareSize*(boardSize(1)-2)/2; % y-offset of frame o wrt frame f
H_o2f = Tx(x_o2f)*Ty(y_o2f)*Rx(pi);

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
% [ALLOWS RUNNING THIS SECTION ONLY] vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
load('Lab07_calibrationParameters.mat');
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

% [ALLOWS RUNNING THIS SECTION ONLY]
% Open Wall-E visualization and get figure handle
if ~exist('figWallE','var') || ~ishandle(figWallE)
    figWallE = open('Lab07_Wall-E.fig');
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

% Recover the patch object from the Wall-E visualization figure handle
p_b = ptc_b.Vertices.'; % <-- Note the transpose
% Do not make the points homogeous, this will make projectWithFalseDepth
% run a little faster.

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
% [ALLOWS RUNNING THIS SECTION ONLY] vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
load('Lab07_calibrationParameters.mat');
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

% [ALLOWS RUNNING THIS SECTION ONLY]
% Open Wall-E visualization and get figure handle
if ~exist('figWallE','var') || ~ishandle(figWallE)
    figWallE = open('Lab07_Wall-E.fig');
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

% Recover the patch object from the Wall-E visualization figure handle
p_b = ptc_b.Vertices.'; % <-- Note the transpose
% Do not make the points homogeous, this will make projectWithFalseDepth
% run a little faster.

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

% /////////////////////////////////////////////////////////////////////////
% NEW CODE
% Initialize an MPEG-4 video writer object
vid = VideoWriter('Lab07_WallE_Driving.mp4','MPEG-4');
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
    
    %{
    % /////////////////////////////////////////////////////////////////////
    % NEW CODE (figure handle)
    % "Take a picture" of the figure handle
    % (this assumes "fig" is the figure object you want in the video)
    frame = getframe(fig);
    % Write the frame to the video
    writeVideo(vid,frame);
    % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    %}
    
    % /////////////////////////////////////////////////////////////////////
    % "Take a picture" of the axes handle
    % (this assumes "axs" is the axes object you want in the video)
    frame = getframe(axs);
    % Write the frame to the video
    writeVideo(vid,frame);
    % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
end

% /////////////////////////////////////////////////////////////////////////
% NEW CODE
% Close video writer object
close(vid);
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

%% Real-time example
% [ALLOWS RUNNING THIS SECTION ONLY] vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
load('Lab07_calibrationParameters.mat');
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

% [ALLOWS RUNNING THIS SECTION ONLY]
% Open Wall-E visualization and get figure handle
if ~exist('figWallE','var') || ~ishandle(figWallE)
    figWallE = open('Lab07_Wall-E.fig');
end

% Recover the patch object from the Wall-E visualization figure handle
ptc_b = findobj(figWallE,'Type','patch','Tag','Wall-E');

% Recover the patch object from the Wall-E visualization figure handle
p_b = ptc_b.Vertices.'; % <-- Note the transpose
% Do not make the points homogeous, this will make projectWithFalseDepth
% run a little faster.

% Define checkerboard info
[boardSize,squareSize] = checkerboardPoints2boardSize( cameraParams.WorldPoints );
% Define Frame o relative to the fiducial frame f
x_o2f = squareSize*(boardSize(2)-2)/2; % x-offset of frame o wrt frame f
y_o2f = squareSize*(boardSize(1)-2)/2; % y-offset of frame o wrt frame f
H_o2f = Tx(x_o2f)*Ty(y_o2f)*Rx(pi);

% Define Wall-E’s body-fixed frame relative to frame o. 
% Note: We will initially use the identity, but we can change this frame
%       if we want Wall-E to drive.
H_b2o = eye(4);

% Create a figure for the real-time example
im = get(prv,'CData');              % Get a new image from the camera
figRealTime = figure('Name','"Real-time Example"'); % Create a new figure
axsRealTime = axes('Parent',figRealTime);      % Create an axes in the figure
imgRealTime = imshow(im,'Parent',axsRealTime); % Show the image in the axes
hold(axsRealTime,'on');
addSingleLight(axsRealTime);
ptc_m = copyobj(ptc_b,axsRealTime);

% Define square size for checkerboard fiducial
squareSize = 19.05; % <-- Confirm that this is correct for your checkerboard

% Create a loop to recover extrinsics
while ishandle(figRealTime)
    % Allow preview and figure(s) to update
    drawnow
    
    % Get image from preview
    im = get(prv,'CData');

    % Define p_m from image of checkerboard
    % NOTE: MATLAB's definition of "imagePoints" relates to p_m as follows:
    %       % Define "imagePoints" from p_m
    %       imagePoints = p_m(1:2,:).'; % <-- Note the transpose
    %       % Define p_m from "imagePoints"
    %       p_m(1:2,:) = imagePoints.'; % <-- Note the transpose
    %       p_m(3,:) = 1;               % <-- Convert to homogeneous, 2D 
    %                                         pixel position
    [imagePoints,boardSize] = detectCheckerboardPoints(im);
    
    % Check if full checkerboard was detected
    if any(~isfinite(imagePoints),'all') || any(boardSize == 0)
        % One or more checkerboard points is not tracked
        % Hide Wall-E
        set(ptc_m,'Visible','off');
        % Continue to next iteration of while-loop
        continue
    else
        % Checkerboard is tracked
        % Show Wall-E
        set(ptc_m,'Visible','on');
    end
    
    % Define p_f given "boardSize" and "squareSize"
    % NOTE: MATLAB's definition of "worldPoints" relates to p_f as follows:
    %       % Define "worldPoints" from p_f
    %       worldPoints = p_f(1:2,:).'; % <-- Note the transpose
    %       % Define p_f from "worldPoints"
    %       p_f(1:2,:) = worldPoints.'; % <-- Note the transpose
    %       p_f(3,:) = 0;               % <-- Define z-coordinate
    %       p_f(4,:) = 1;               % <-- Convert to homogeneous, 3D
    %                                   %     coordinate relative to the
    %                                   %     fiducial frame
    [worldPoints] = generateCheckerboardPoints(boardSize,squareSize);
    
    % Recover the checkerboard pose relative to the camera frame (H_f2c)
    [R_c2f, tpose_d_f2c] = extrinsics(...
        imagePoints,worldPoints,cameraParams);
    R_f2c = R_c2f.'; 
    d_f2c = tpose_d_f2c.';
    H_f2c = [R_f2c, d_f2c; 0,0,0,1];
    
    % Define applicable extrinsics
    H_b2c = H_f2c*H_o2f*H_b2o;
    
    % Define projection
    P_b2m = A_c2m * H_b2c(1:3,:);
    
    % Project points with added false depth
    p_m_falseDepth = projectWithFalseDepth(p_b,P_b2m);

    % Update Image 
    set(imgRealTime,'CData',im);
    
    % Update Wall-E
    ptc_m.Vertices = p_m_falseDepth.';
    drawnow;
end


