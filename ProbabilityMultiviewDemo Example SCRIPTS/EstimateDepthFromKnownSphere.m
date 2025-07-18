%% EstimateDepthFromKnownSphere
% Estimate the distance of a sphere of known radius given a camera with
% known intrinsics
%
%   NOTE: This function uses camera parameters for MATLAB 2023b or older. 
%         If you are using a newer version of MAT
%   M. Kutzer, 23Mar2021, USNA

clear all
close all
clc

%% Define makeVIDEO
makeVIDEO = false;

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

%% Define intrinsic matrix
A_c2m = cameraParams.IntrinsicMatrix.';

%% Create figures
% -> 3D Figure
% Initialize figure and axes
fig3D = figure('Name','3D View','Color',[1 1 1]);
axs3D = axes('Parent',fig3D,'NextPlot','Add','DataAspectRatio',[1 1 1]);
addSingleLight(axs3D);
view(axs3D,3);
% Plot camera
plt_c = plotCamera('Size',50,'Color','b');

% -> 2D Figure
% Define image size
x_m_MAX = cameraParams.ImageSize(2);
y_m_MAX = cameraParams.ImageSize(1);
% Initialize figure and axes
fig2D = figure('Name','Image Generator','Units','Pixels',...
    'Position',[10,10,x_m_MAX,y_m_MAX],'Color',[1 1 1]);
centerfig(fig2D);
axs2D = axes('Parent',fig2D,'Units','Normalized','Position',[0,0,1,1],...
    'Visible','off','NextPlot','Add','DataAspectRatio',[1 1 1]);
xlim(axs2D,0.5 + [0,x_m_MAX]);
ylim(axs2D,0.5 + [0,y_m_MAX]);
set(fig2D,'Visible','off');


%% Define sphere
r = 51/2;   % 51 mm diameter sphere
sfit.Radius = r;
sfit.Center = [0,0,0];

%% Create patch
p = patchSphere(sfit,2000);
ptc3D = patch(p,'Parent',axs3D,'FaceColor','r','EdgeColor','None');
ptc2D = patch(p,'Parent',axs2D,'FaceColor','r','EdgeColor','None');

%% Define rigid body transform to move 3D sphere
h_s2c = triad('Parent',axs3D,'Scale',55,'AxisLabels',{'x_s','y_s','z_s'});
set(ptc3D,'Parent',h_s2c);

%% Show Simulated Image
% Get image simulated image
drawnow
frm = getframe(axs2D);
im = frm.cdata;
im = imresize(im,[y_m_MAX,x_m_MAX]);
% Visualize simulated image
figIM = figure('Name','Simulated Image');
imgIM = imshow(im);
axsIM = get(imgIM,'Parent');
set(axsIM,'NextPlot','add','Visible','On');
xlabel(axsIM,'x (pixels)');
ylabel(axsIM,'y (pixels)');
% Define segmentation and centroid overlay
ptc_bm(1) = patch('Vertices',nan(3,2),'Faces',1:3,'Parent',axsIM,...
    'EdgeColor','c','FaceColor','none','FaceAlpha',0.1,'LineWidth',2.0,...
    'LineStyle','-');
ptc_bm(2) = patch('Vertices',nan(3,2),'Faces',1:3,'Parent',axsIM,...
    'EdgeColor','k','FaceColor','k','FaceAlpha',0.1,'LineWidth',1.2,...
    'LineStyle',':');
plt_cm = plot(axsIM,nan,nan,'xc','MarkerSize',10,'LineWidth',1.2);

%% Test
z_c_MAX = 10.0*12*25.4;   % 18.0 feet, converted to mm
z_c_MIN =  0.5*12*25.4;   %  0.5 feet, converted to mm
n = 200;
z_c_TRU = linspace(z_c_MAX,z_c_MIN,n);

% Define vertices in the sphere frame
v_s = p.Vertices.';
v_s(4,:) = 1;

% Move 3D ball
set(h_s2c,'Matrix',Tz(z_c_MAX));
xlim(axs3D,[-200,200]);
ylim(axs3D,[-200,200]);
zlim(axs3D,[-1.5*r,z_c_MAX + 1.5*r]);

if makeVIDEO
    vid2D = VideoWriter('Ball, 2D View.mp4','MPEG-4');
    vid3D = VideoWriter('Ball, 3D View.mp4','MPEG-4');
    open(vid2D);
    open(vid3D);
end

iter = 0;
for z_c = z_c_TRU
    iter = iter+1;
    % Define extrinsics
    H_s2c = Tz(z_c);
    
    % Transform patches
    % -> 3D
    v_c = H_s2c * v_s;
    set(h_s2c,'Matrix',H_s2c);
    % -> 2D
    sv_m = A_c2m*v_c(1:3,:);
    z_cs = sv_m(3,:);
    v_m = sv_m./z_cs;
    v_m(3,:) = z_cs;
    set(ptc2D,'Vertices',v_m.');
    
    % Stop if ball leaves FOV
    if ...
            ( nnz(v_m(1,:) < 0) || nnz(v_m(1,:) > y_m_MAX) ) ||...
            ( nnz(v_m(2,:) < 0) || nnz(v_m(2,:) > x_m_MAX) )
        break
    end
    
    % Get image simulated image
    drawnow
    frm = getframe(axs2D);
    im = frm.cdata;
    im = imresize(im,[y_m_MAX,x_m_MAX]);
    
    % Update image simulation
    set(imgIM,'CData',im);

    % Segment ball from image
    bin = segmentRedBall(im);
    a_EST(1,iter) = bwarea(bin);
    [y_m,x_m] = bwCentroid(bin);

    % Visualize segmentation and centroid
    if nnz(bin) > 0
        set(plt_cm,'XData',x_m,'YData',y_m,'Visible','on');
        % -> Trace the boundaries
        bnd_m = bwboundaries(bin);
        set(ptc_bm,'Vertices',fliplr(bnd_m{1}),'Faces',1:size(bnd_m{1},1),'Visible','on');
    else
        set(plt_cm,'Visible','off');
        set(ptc_bm,'Visible','off');
    end

    % Recover depth and estimate x/y position
    z_c_EST(1,iter) = depthFromSphereArea(A_c2m,r,a_EST(1,iter));
    a_TRU(1,iter) = sphereAreaFromDepth(A_c2m,r,z_c);
    X_m_EST(:,iter) = [y_m; x_m; 1];
    X_c_EST(:,iter) = z_c_EST(1,iter) * (A_c2m^(-1)) * X_m_EST(:,iter);
    % Recover "true" center
    X_c_TRU(:,iter) = H_s2c(1:3,4);
    
    if makeVIDEO
        frm = getframe(fig3D);
        writeVideo(vid3D,frm);
        frm = getframe(fig2D);
        writeVideo(vid2D,frm);
    end
end
if makeVIDEO
    close(vid2D);
    close(vid3D);
end
z_c_TRU(iter:end) = [];

%% PLOT RESULTS
figure;
axes;
hold on
plot(z_c_TRU,z_c_TRU - z_c_EST);
xlabel('True z_c (mm)');
ylabel('(True z_c) - (Estimated z_c) (mm)');
if makeVIDEO
    drawnow
    saveas(gcf,'z Error vs z_c.png','png');
end

figure;
axes;
hold on
err = sqrt( sum( (X_c_TRU - X_c_EST).^2, 1 ) );
plot(z_c_TRU,err);
xlabel('True z_c (mm)');
ylabel('norm(X^{c}_{TRU} - X^{c}_{EST}) (mm)');
if makeVIDEO
    drawnow
    saveas(gcf,'Total Error vs z_c.png','png');
end

figure;
axes;
hold on
plot(a_EST,z_c_TRU,'g');
plot(a_EST,z_c_EST,'m');
xlabel('Area (pixels)');
ylabel('z_c (mm)');
if makeVIDEO
    drawnow
    saveas(gcf,'Area vs z_c.png','png');
end

figure;
axes;
hold on
plot(a_EST,a_TRU);
xlabel('Estimated Area From Image (pixels)');
ylabel('Area from Z_c (pixels)');
if makeVIDEO
    drawnow
    saveas(gcf,'Estimated area vs actual area','png');
end

%% Fit area correction
pAest2Atru = polyfit(a_EST,a_TRU,1);
save('AreaCorrection.mat','pAest2Atru');