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

%% Visualize intrinsic uncertainty
fig3D = figure('Name','Intrinsic Uncertainty');
axs3D = axes('Parent',fig3D,'NextPlot','add');
view(axs3D,3);
xlabel(axs3D,'x (pixels)');
ylabel(axs3D,'y (pixels)');
zlabel(axs3D,'Image Number');

% Define unit square in x/y camera frame
p_c = [...
    -1,1,  1,-1;...
    -1,-1, 1, 1;...
     1, 1, 1, 1]*0.5;
for i = 1:nImages
    % Skip empty values
    if isempty(A_c2m_i{i})
        continue
    end

    % Project points
    tilde_p_m = A_c2m*p_c;
    p_m = tilde_p_m./tilde_p_m(3,:);

    % Offset "z" using image index (for visualization only)
    p_m_i = p_m;
    p_m_i(3,:) = i;
    
    % Calculate likelihood associated with sample
    y = mvnpdf(Av_c2m_i(:,i).',barAv_c2m.',sigAv_c2m);

    % Visualize intrinsic sample
    ptc3D_i(i) = patch('Parent',axs3D,'Vertices',p_m_i.','Faces',1:4,...
        'EdgeColor','b','FaceColor','b','FaceAlpha',0.2);
    plt3D_i(i) = plot3(axs3D,A_c2m_i{i}(1,3),A_c2m_i{i}(2,3),i,'*m');
    txt3D_i(i) = text(axs3D,A_c2m_i{i}(1,3),A_c2m_i{i}(2,3),i,...
        sprintf('%6.4f',y));
end