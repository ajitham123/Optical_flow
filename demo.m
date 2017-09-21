%% Optical flow and FOE for artificial images
%   clear all; close all;
im_scale = false;
if_sub_pixel = false;
real_FOE = [200,400];

% load images:
I1 = imread('images/office.jpg');
% I2 = imread('images/indoor_foe_180_315.jpg');
I2 = FOE_db(I1,real_FOE,1.015);
if im_scale
    I1 = imresize(I1,0.4);
    I2 = imresize(I2,0.4);
end
FOE = [0; 0];
for i=1:1
% determine the individual flow vectors:
graphics = true;
[point, flow_mag, angle] = get_optical_flow_edges(I1, I2, graphics,if_sub_pixel);

% find FOE
 FOE = FOE+find_FOE(point, flow_mag, angle,graphics);
end
FOE=FOE/1
norm(real_FOE'-FOE)

%% Optical flow and FOE for a video

% clc; clear all; close all;
% graphics = true;
% if_sub_pixel = false;
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
%     [point, flow_mag, angle] = get_optical_flow_edges(I1, I2, graphics,if_sub_pixel);
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
if_sub_pixel = false;
graphics = false;
stride = 50;

% load images:
I1 = imread('images/street.jpg');
% I2 = imread('images/indoor_foe_180_315.jpg');
h_map = zeros(floor(size(I1,1)/stride),floor(size(I1,2)/stride));

for i=stride:stride:size(I1,1)
    for j=stride:stride:size(I1,2)
        real_FOE = [i,j];
        I2 = FOE_db(I1,real_FOE,1.015);
        
        % determine the individual flow vectors:
        [point, flow_mag, angle] = get_optical_flow_edges(I1, I2, graphics,if_sub_pixel);
        
        % find FOE
        FOE = find_FOE(point, flow_mag, angle,graphics);
        
        h_map(i/stride,j/stride) = norm(real_FOE'-FOE);
        fprintf("\n i = %d, j = %d, FOE error = %f",i,j,h_map(i/stride,j/stride));
    end
end
fprintf("\n Mean error = %f ",mean(h_map(:)));
heatmap(h_map)


