%% ExampleIntrinsicUncertainty
% This script estimates intrinsic uncertainty given calibration images of a
% known checkerboard.
%
%   NOTE: This script relies on data saved from ExampleCameraCalibration.m
%
%   M. Kutzer, 18Jul2024, USNA

clear all
close all
clc

%% Load saved data
load('ExampleCalibrationParameters.mat');

%% Define checkerboard info
[boardSize,squareSize] = checkerboardPoints2boardSize( cameraParams.WorldPoints );

%% Parse mean intrinsics
barA_c2m = cameraParams.IntrinsicMatrix.'; % <-- Note the transpose
% Define vector form
barAv_c2m = veeIntrinsics(barA_c2m);

%% Calculate intrinsic samples
A_c2m_i = {};
Av_c2m_i = [];
for i = 1:nImages
    % Define filename
    imName = sprintf('%s%03d.png',imBaseName,imagesUsed(i));
    % Load image
    im = imread( fullfile(calFolderName,imName) );

    % Calculate intrinsics sample
    A_c2m_i{i} = imageToIntrinsics(im,cameraParams,squareSize,boardSize);
    % Define vector form
    if ~isempty(A_c2m_i{i})
        Av_c2m_i(:,i) = veeIntrinsics(A_c2m_i{i});
    else
        Av_c2m_i(:,i) = nan(5,1);
    end
end

% Flag empty elements
tfIsEmpty = cellfun(@isempty,A_c2m_i);

%% Define intrinsic covariance
sigAv_c2m = ...
    covGivenMean(Av_c2m_i(:,~tfIsEmpty).',barAv_c2m.'); % <-- Note transpose

%% Visualize intrinsic samples

% Initialize figure and axes
fig3D = figure('Name','Intrinsic Samples');
axs3D = axes('Parent',fig3D,'NextPlot','add');
view(axs3D,3);
xlabel(axs3D,'x^c / z^c (unitless)');
ylabel(axs3D,'y^c / z^c (unitless)');
zlabel(axs3D,'Image Number');

% Define image x/y limits in pixels
xLims_m = [0,cameraParams.ImageSize(2)];
yLims_m = [0,cameraParams.ImageSize(1)];
% Define image bounds
p_m = [...
    xLims_m([1,2,2,1]);...
    yLims_m([1,1,2,2])];
p_m(3,:) = 1;

% Visualize image bounds as scaled camera coordinates
for i = 1:nImages
    % Skip empty values
    if isempty(A_c2m_i{i})
        continue
    end

    % Define scaled camera-referenced points
    tilde_p_c = ( A_c2m_i{i} )^(-1)*p_m;

    % Offset "z" using image index (for visualization only)
    tilde_p_c(3,:) = i;
    
    % Calculate likelihood associated with sample
    y = mvnpdf(Av_c2m_i(:,i).',barAv_c2m.',sigAv_c2m);

    % Visualize intrinsic sample
    ptc3D_i(i) = patch('Parent',axs3D,'Vertices',tilde_p_c.','Faces',1:4,...
        'EdgeColor',rand(1,3),'FaceColor','none','FaceAlpha',0.2);
    txt3D_i(i) = text(axs3D,...
        mean(tilde_p_c(1,:)),...
        mean(tilde_p_c(2,:)),...
        mean(tilde_p_c(3,:)),...
        sprintf('%6.5f',y));
end

%% Calculate and visualize reprojection errors
barErr = [];
smpErr = [];
for i = 1:nImages
    % Define filename
    imName = sprintf('%s%03d.png',imBaseName,imagesUsed(i));
    % Load image
    im = imread( fullfile(calFolderName,imName) );

    % Error with mean intrinsics
    barErr(i) = imageToReprojectionError(im,cameraParams,squareSize,barA_c2m);

    % Skip empty values
    if isempty(A_c2m_i{i})
        smpErr(i) = nan;
        continue
    end

    % Error with sample intrinsics
    smpErr(i) = imageToReprojectionError(im,cameraParams,squareSize,A_c2m_i{i});
end

% Plot results
fig = figure('Name','Reprojection Error Comparison');
axs = axes('Parent',fig,'NextPlot','add');
xlabel(axs,'Image Number');
ylabel(axs,'RMS Reprojection Error (pixels)');
plt_b = bar(barErr,'Parent',axs);
plt_s = bar(smpErr,'Parent',axs);
legend(axs,'Calibration Intrinsics','Sample Intrinsics');

%% Visualize intrinsic uncertainty
% Define image x/y limits in pixels
xLims_m = [0,cameraParams.ImageSize(2)];
yLims_m = [0,cameraParams.ImageSize(1)];
% Define meshgrid sampling of the image space
d = 10;
m = round( diff(xLims_m)/d ); % x-samples
n = round( diff(yLims_m)/d ); % y-samples
[X_m,Y_m] = meshgrid(...
    linspace(xLims_m(1),xLims_m(2),m),...
    linspace(yLims_m(1),yLims_m(2),n) );
% Reshape to define pixel coordinates 
p_m = [];
p_m(1,:) = reshape(X_m,1,[]);
p_m(2,:) = reshape(Y_m,1,[]);
p_m(3,:) = 1;

% Generate random intrinsic matrix sampling
k = 1000; % intrinsic samples
Av_c2m_k = mvnrnd(barAv_c2m.',sigAv_c2m,k).';

% Calculate mean scaled camera coordinates
tilde_pBar_c = ( barA_c2m )^(-1)*p_m;

% Calculate scaled camera coordinates
tilde_dp_c_k = zeros(k,m*n);
for i = 1:k
    % Define scaled camera-referenced points
    tilde_p_c = ( wedgeIntrinsics(Av_c2m_k(:,i)) )^(-1)*p_m;

    % Define zero-mean points
    tilde_dp_c = tilde_p_c - tilde_pBar_c;

    % Define distance to mean
    tilde_dp_c_k(i,:) = sqrt( sum( tilde_dp_c.^2,1 ) );
end

% Define variance in scaled camera coordinates
tilde_dpVar_c = var(tilde_dp_c_k);

% Plot results
fig = figure('Name','Intrinsic Uncertainty');
z_max = max(tilde_dpVar_c);
axs = axes('Parent',fig,'NextPlot','add','YDir','reverse');%,...
    %'DataAspectRatio',(1/z_max)*[1 1 z_max]);
xlim(axs,xLims_m);
ylim(axs,yLims_m);
zlim(axs,[0,max(tilde_dpVar_c)]);
xlabel(axs,'x^m (pixels)');
ylabel(axs,'y^m (pixels)');

% 3D points (for debugging)
%plt = plot3(axs,p_m(1,:),p_m(2,:),tilde_dpVar_c,'.');

% Surface Plot
Z = reshape(tilde_dpVar_c,n,m);
srf = surf(axs,X_m,Y_m,Z,'EdgeColor','none');
cbr = colorbar(axs,'eastoutside');

%{
% Add z & color labels
% -> label string
str = '$$\sigma^2\left( \left| \frac{p^c}{z^c} - \frac{\bar{p}^c}{\bar{z}^c} \right|\right)$$';
% -> z-label
zlabel(axs,str,'Interpreter','Latex');
% -> color label
cbr.Label.Interpreter = 'Latex';
cbr.Label.String = str;
%}
