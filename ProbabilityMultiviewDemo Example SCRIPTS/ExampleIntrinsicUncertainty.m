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
% Vector form
barAv_c2m = vee(barA_c2m);

%% Calculate intrinsic samples
