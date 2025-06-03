%% DepthError_2Camera_ChangePose
% Calculate the depth error associated with error in individual point
% locations. Evaluate the change in the error relationship as the cameras
% are moved farther appart while keeping a common z_c between cameras.
%
%   NOTE: This script uses several simple internal functions to improve
%         readability. Scroll to the bottom of the script to review these
%         functions as needed.
%
%   M. Kutzer, 28Mar2022, USNA

clear all
close all
clc

makeVideo = false;

if makeVideo
    % Setup to create video
    vidTitle = sprintf('DepthError_2Camera_3D.mp4');
    vidObj_3D = VideoWriter(vidTitle,'MPEG-4');
    open(vidObj_3D);

    % Setup to create video
    vidTitle = sprintf('DepthError_2Camera_cA.mp4');
    vidObj_c(1) = VideoWriter(vidTitle,'MPEG-4');
    open(vidObj_c(1));

    % Setup to create video
    vidTitle = sprintf('DepthError_2Camera_c2.mp4');
    vidObj_c(2) = VideoWriter(vidTitle,'MPEG-4');
    open(vidObj_c(2));
end

%% Create figure & axes
fig = figure('Name','Estimate Depth Error');
axs = axes('Parent',fig);
hold(axs,'on');
daspect(axs,[1 1 1]);
view(axs,3);
xlabel(axs,'x');
ylabel(axs,'y');
zlabel(axs,'z');
addSingleLight(axs);

% Adjust axes limits
xlim(axs,[-600, 600]);
ylim(axs,[-900, 900])
zlim(axs,[-200,1200]);

% Figure for zoomed error
figZoom = figure('Name','Zoomed Error');
axsZoom = axes('Parent',figZoom);
hold(axsZoom,'on');
daspect(axsZoom,[1 1 1]);
view(axsZoom,3);
xlabel(axsZoom,'x (mm)');
ylabel(axsZoom,'y (mm)');
zlabel(axsZoom,'z (mm)');
addSingleLight(axsZoom);

% Set figure position and color for video
set(fig,'Color',[1 1 1],...
    'Units','Inches','Position',[0.25,0.65,6.27,4.70]);
set(figZoom,'Color',[1 1 1],...
    'Units','Inches','Position',[6.80,0.65,6.27,4.70]);

%% Define camera intrinsics
% NOTE: For this error approximation, we will *not* use a full camera
%       model. The cameras are assumed to be "perfect pinhole cameras."
%       Introducing lens effects etc. will further increase error.
A_ca2ma = [...
    460, 0, 320;...
    0, 470, 180;...
    0,   0,   1];

A_cb2mb = [...
    470, 0, 332;...
    0, 465, 176;...
    0,   0,   1];

%% Initialize plots of cameras and FOVs
h_ca2a = drawDFKCam(axs,eye(4));
h_cb2a = drawDFKCam(axs,eye(4));

res = [640,360];
fov_c(1) = plotCameraFOV(h_ca2a,A_ca2ma,1300,res);
set(fov_c(1),'FaceAlpha',0.1,'FaceColor','r');

fov_c(2) = plotCameraFOV(h_cb2a,A_cb2mb,1300,res);
set(fov_c(2),'FaceAlpha',0.1,'FaceColor','b');

%% Create image visualizations
for i = 1:2
    fig2D(i) = figure('Name',sprintf('Simulated Image %d',i));
    axs2D(i) = axes('Parent',fig2D(i));
    hold(axs2D(i),'on');
    daspect(axs2D(i),[1 1 1]);
    xlim(axs2D(i),[0,res(1)]+[-0.5,0.5]);
    ylim(axs2D(i),[0,res(2)]+[-0.5,0.5]);
    ttl(i) = title(axs2D(i),sprintf('Camera %d',i));
    plt2D_tru(i) = plot(axs2D(i),nan,nan,'xk');
    plt2D_rec(i) = plot(axs2D(i),nan,nan,'+k');
    plt2D_com(i) = plot(axs2D(i),nan,nan,'ok');
    lgnd2D(i) = legend(axs2D(i),[plt2D_tru(i), plt2D_rec(i),plt2D_com(i)],...
        'True Center',...
        'Recovered Center',...
        'Common Center');
end
set(axs2D,'ydir','reverse');
xlabel(axs2D,'x (pixels)');
ylabel(axs2D,'y (pixels)');

%% Initialize error region plot
ptc_err = patch('Vertices',nan(3,3),'Faces',1:3,'Parent',axs,...
    'EdgeColor','none','FaceColor','b','FaceAlpha',0.6);
ptc_err_z = copyobj(ptc_err,axsZoom);

%% Define "true" point
p_a_true = [0;0;0;1];
plt_a = plot3(axs,p_a_true(1),p_a_true(2),p_a_true(3),'*m',...
    'MarkerSize',10,'LineWidth',1);
plt_a_z = copyobj(plt_a,axsZoom);

%% Define test parameters
% Distance from point to each camera (true z_c for both)
s = 1100;
% Angle of seperation between cameras
k = 60;
alpha = deg2rad( linspace(0,180,k+2) );
alpha(1) = [];
alpha(end) = [];

% Error magnitude (i.e. radius)
n = 15; % Number of error radius samples
r = linspace(0,10,n);       % pixels

% Error angle
m = 15; % Number of sampling angles
phi = linspace(0,2*pi,m+1); % radians
phi(end) = [];

a = 1;
%for a = 1:numel(alpha)
%tic;
% -> Define camera locations & projection matrices
H_ca2a = Rx( alpha(a)/2)*Tz(s)*Ry(pi);
H_cb2a = Rx(-alpha(a)/2)*Tz(s)*Ry(pi);

% Define extrinsics
H_a2ca = invSE(H_ca2a);
H_a2cb = invSE(H_cb2a);

% Define projection matrices
P_a2ma = A_ca2ma*H_a2ca(1:3,:);
P_a2mb = A_cb2mb*H_a2cb(1:3,:);

% Define relative extrinsics
H_cb2ca = H_a2ca*H_cb2a;

% Update 3D view
set(h_ca2a,'Matrix',H_ca2a);
set(h_cb2a,'Matrix',H_cb2a);

%% Project points defining "true" pixel coordinates
[p_ma,z_ca_true] = projectPoints(p_a_true, P_a2ma);
[p_mb,z_cb_true] = projectPoints(p_a_true, P_a2mb);

%% Define radius error & sampling angles
p_a_est = nan(3,n,n,m,m);
ssd_err = nan(n,n,m,m);
for ia = 1:n
    for ib = 1:n
        for ja = 1:m
            for jb = 1:m
                % Define error vector
                err_ma = r(ia)*[cos(phi(ja)); sin(phi(ja)); 0];
                err_mb = r(ib)*[cos(phi(jb)); sin(phi(jb)); 0];

                % Recover depth
                [z_ca,z_cb] = recoverDepth(...
                    p_ma+err_ma,p_mb+err_mb,A_ca2ma,A_cb2mb,H_cb2ca);

                % Recover point in 3D
                p_a1 = recoverPoints(p_ma+err_ma, z_ca, P_a2ma);
                p_a2 = recoverPoints(p_mb+err_mb, z_cb, P_a2mb);

                % Average result
                p_a = mean([p_a1,p_a2],2);
                
                % Project common point back into the image
                p_ma_com = projectPoints(p_a, P_a2ma);
                p_mb_com = projectPoints(p_a, P_a2mb);
                % Append p_a
                %p_a_est(:,ia,ib,ja,jb) = p_a(1:3,:);

                % Calculate error
                %ssd_err(ia,ib,ja,jb) = norm(p_a - p_a_true);

                % Update plots
                set(plt2D_tru(1),'XData',p_ma(1),'YData',p_ma(2));
                set(plt2D_tru(2),'XData',p_mb(1),'YData',p_mb(2));
                set(plt2D_rec(1),...
                    'XData',p_ma(1)+err_ma(1),...
                    'YData',p_ma(2)+err_ma(2));
                set(plt2D_rec(2),...
                    'XData',p_mb(1)+err_mb(1),...
                    'YData',p_mb(2)+err_mb(2));
                set(plt2D_com(1),'XData',p_ma_com(1),'YData',p_ma_com(2));
                set(plt2D_com(2),'XData',p_mb_com(1),'YData',p_mb_com(2));

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
        end
    end
end

%% Plot 3D error & statistics
%p_a_err = p_a_est - repmat(p_a_true(1:3,:),1,n,n,m,m);
%p_a_err = p_a_est - repmat(p_a_true(1:3,:),1,m*m,n,n);
%rshp_p_a_est = reshape(p_a_est,3,[]);
%faces = convhull(rshp_p_a_est(1,:),rshp_p_a_est(2,:),rshp_p_a_est(3,:));
%set(ptc_err,'Vertices',rshp_p_a_est.','Faces',faces);
%set(ptc_err_z,'Vertices',rshp_p_a_est.','Faces',faces);
%drawnow


%fname = sprintf('ErrData%.2f,%d,%d,%d,%02d.mat',s,k,n,m,a);
%save(fname,'p_a_est','ssd_err','s','k','m','n','alpha','r','phi','a')
%toc
%end

if makeVideo
    % Close video obj
    close(vidObj_3D);
    close(vidObj_c(1));
    close(vidObj_c(2));
end


%% Define internal utility functions
function [p_m,z_c] = projectPoints(p_a, P_a2m)
% Calculate scaled pixel coordinate(s)
tilde_p_m = P_a2m*p_a;

% Recover depth(s)
z_c = tilde_p_m(3,:);

% Recover homogeneous pixel coordinate
p_m = tilde_p_m./z_c;
end

function p_a = recoverPoints(p_m, z_c, P_a2m)
% Isolate "rotation" and "translation" portions of projection matrix
Q_a2m = P_a2m(1:3,1:3);
V_a2m = P_a2m(1:3,4);

% Recover 3D point coordinate(s) & make homogeneous
p_a = (Q_a2m^(-1))*(z_c.*p_m - V_a2m);
p_a(4,:) = 1;
end

function [z_ca,z_cb] = recoverDepth(p_ma,p_mb,A_ca2ma,A_cb2mb,H_cb2ca)
% Define scaled point coordinates relative to camera frame
tilde_p_ca = (A_ca2ma^(-1))*p_ma;
tilde_p_cb = (A_cb2mb^(-1))*p_mb;

% Isolate rotation & translation of extrinsics
R_cb2ca = H_cb2ca(1:3,1:3);
d_cb2ca = H_cb2ca(1:3,4);

% Calcualte depth
z_c_ALL = pinv( [tilde_p_ca, -R_cb2ca*tilde_p_cb] )*d_cb2ca;

% Isolate depth parameters
z_ca = z_c_ALL(1,:);
z_cb = z_c_ALL(2,:);
end
