%% VisualizeRigidBodyTransform
% This script provides a simple example of defining the postion and 
% orientation of two frames relative to one another using motion 
% primitives. This script also demonstrates the visualization of two frames 
% using the "triad" function. 
%
%   M. Kutzer, 30Jul2020, USNA

%% Create figure & axes
fig = figure('Name','Example Visualize Rigid Body Transformation');
axs = axes('Parent',fig);
daspect(axs,[1 1 1]);
hold(axs,'on');

%% Adjust the axes properties to pretty-up the plot
axis(axs,[-10,110,-10,50,-10,10]);
grid(axs,'on');
set(axs,'XTick',[-10:10:110]);

%% Create coordinate frame representations
h_1to0 = triad('Parent',axs,'AxisLabels',{'x_1','y_1','z_1'},...
    'Scale',12,'LineWidth',2);
h_2to1 = triad('Parent',h_1to0,'AxisLabels',{'x_2','y_2','z_2'},...
    'Scale',12,'LineWidth',2);

%% Update pose of Frame 2 relative to Frame 1
H_2to1 = Tx(96)*Ty(25)*Ry(-pi/2);   % Define mathematically
set(h_2to1,'Matrix',H_2to1);        % Update visualization
