function a = sphereAreaFromDepth(A_c2m,r,z_c)
% SPHEREAREAFROMDEPTH calculates the pixel area of a sphere of known 
% radius given camera extrinsics and depth (z_c) from the camera.
%   z_c = SPHEREAREAFROMDEPTH(A_c2m,r,z_c)
%
%   Input(s) 
%       A_c2m - 3x3 camera intrinsics
%       r     - scalar radius of sphere (units must calibration)
%       z_c   - Depth (or distance) from camera relative to the camera
%               frame
%
%   Output(s)
%       a     - scalar area of the sphere in the image (pixels)
%
%   Note: This function currently assumes an intrinsic matrix with zero 
%         shear (i.e. having the following form):
%                 [sx,  0, u]
%         A_c2m = [ 0, sy, v]
%                 [ 0,  0, 1]
%       
%   M. Kutzer, 07Apr2021, USNA

%% Check input(s)
% TODO - check inputs

%% Calculate z_c
sx = A_c2m(1,1);
sy = A_c2m(2,2);
a = ( (pi*abs(r*sx)*abs(r*sy))./(z_c.^2) );