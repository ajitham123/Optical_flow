%% Optical flow and FOE for 2 images
% clc; clear all; close all;
im_scale = false;

% load images:
I1 = imread('images/shapes1_zoom.png');
I2 = imread('images/shapes1_zoom1.png');
if im_scale
    I1 = imresize(I1,.4);
    I2 = imresize(I2,0.4);
end

% determine the individual flow vectors:
graphics = true;
[point, flow_mag, angle] = get_optical_flow_edges(I1, I2, graphics);

% find FOE
FOE = find_FOE(point, flow_mag, angle,graphics);

%% Optical flow and FOE for a video

clc; clear all; close all;
graphics = 1;
im_scale = true;
global frame;

for frame=1:62
    % load images:
    I1 = imread(['C:\Users\Ajith A M\OneDrive\Windows desktop backup\MAVlab\Optical_flow\frames\frames',int2str(frame),'.jpg']);
    I2 = imread(['C:\Users\Ajith A M\OneDrive\Windows desktop backup\MAVlab\Optical_flow\frames\frames',int2str(frame+1),'.jpg']);
    % scale the image if it is too big
    if im_scale
        I1 = imresize(I1,0.4);
        I2 = imresize(I2,0.4);
    end
    
    % determine the individual flow vectors:
    [point, flow_mag, angle] = get_optical_flow_edges(I1, I2, graphics);
%     figure(5);histogram(round(flow_mag(find(flow_mag))))
    
    % find FOE
    FOE = find_FOE(point, flow_mag, angle,graphics);
    pause(.0001);
    frame;
end





