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
barAv_c2m = vee(barA_c2m);

%% Calculate intrinsic samples
A_c2m_i = {};
Av_c2m_i = [];
for i = 1:nImages
    % Define filename
    imName = sprintf('%s%03d.png',imBaseName,imagesUsed(i));
    % Load image
    im = imread( fullfile(calFolderName,imName) );

    % Calculate intrinsics sample
    A_c2m_i{i} = imageToIntrinsics(im,camaraParams,squareSize,boardSize);
    % Define vector form
    if ~isempty(A_c2m_i{i})
        Av_c2m_i(:,i) = vee(A_c2m_i{i});
    else
        Av_c2m_i(:,i) = nan(5,1);
    end
end

% Flag empty elements
tfIsEmpty = cellfun(@isempty,A_c2m_i);

%% Define intrinsic covariance
sigAv_c2m = ...
    covGivenMean(Av_c2m_i(~tfIsEmpty).',barAv_c2m(:,~tfIsEmpty).'); % <-- Note transpose
