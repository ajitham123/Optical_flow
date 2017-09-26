%% Optical flow and FOE for artificial images
%   clear all; close all;
graphics = true;
im_scale = false;
if_sub_pixel = true;
if_ssd = true;
real_FOE = [100,100];

% load images:
I1 = imread('images/street.jpg');
% I2 = imread('images/frames57.jpg');
I2 = FOE_db(I1,real_FOE,1.02);
if im_scale
    I1 = imresize(I1,0.4);
    I2 = imresize(I2,0.4);
end

% determine the individual flow vectors:
[point, flow_mag, angle] = get_optical_flow_edges(I1, I2, graphics,if_sub_pixel,if_ssd);

% find FOE
FOE = find_FOE(point, flow_mag, angle,graphics);

norm(real_FOE'-FOE)

%% Optical flow and FOE for a video

% clc; clear all; close all;
% graphics = true;
% if_sub_pixel = false;
% if_ssd = true;
% im_scale = true;
% global frame;
% 
% for frame=1:62
%     % load images:
%     I1 = imread(['C:\Users\Ajith A M\OneDrive\Windows desktop backup\MAVlab\Optical_flow\frames\frames',int2str(frame),'.jpg']);
%     I2 = imread(['C:\Users\Ajith A M\OneDrive\Windows desktop backup\MAVlab\Optical_flow\frames\frames',int2str(frame+1),'.jpg']);
%     % scale the image if it is too big
%     if im_scale
%         I1 = imresize(I1,0.4);
%         I2 = imresize(I2,0.4);
%     end
%     
%     % determine the individual flow vectors:
%     [point, flow_mag, angle] = get_optical_flow_edges(I1, I2, graphics,if_sub_pixel,if_ssd);
% %     figure(5);histogram(round(flow_mag(find(flow_mag))))
%     
%     % find FOE
%     FOE = find_FOE(point, flow_mag, angle,graphics);
%     pause(.0001);
%     frame;
% end
% 
% 

%% FOE error heatmap of the image 

clc; clear all; close all;
if_sub_pixel = true;
if_ssd = true;
graphics = false;
stride = 50;

% load images:
I1 = imread('images/street.jpg');
% I2 = imread('images/indoor_foe_180_315.jpg');
h_map = zeros(floor(size(I1,1)/stride),floor(size(I1,2)/stride));

for i=stride:stride:size(I1,1)
    for j=stride:stride:size(I1,2)
        real_FOE = [i,j];
        I2 = FOE_db(I1,real_FOE,1.02);
        
        % determine the individual flow vectors:
        [point, flow_mag, angle] = get_optical_flow_edges(I1, I2, graphics,if_sub_pixel,if_ssd);
        
        % find FOE
        FOE = find_FOE(point, flow_mag, angle,graphics);
        
        h_map(i/stride,j/stride) = norm(real_FOE'-FOE);
        fprintf("\n i = %d, j = %d, FOE error = %f",i,j,h_map(i/stride,j/stride));
    end
end
fprintf("\n Mean error = %f ",mean(h_map(:)));
heatmap(h_map)


