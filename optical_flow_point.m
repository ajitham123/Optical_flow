%% Program to find the optical flow along a given direction at a given pixel
%  Author: Ajith Anil Meera, 28th July 2017

%% Evaluate image gradient

clear all; clc; close all;
I1 = rgb2gray(imread('flag1.png'));
I2 = rgb2gray(imread('flag2.png'));

siz = size(I1);
% figure
% imshow(I)
[Gmag, Gdir] = imgradient(I1,'prewitt');

figure
imshowpair(Gmag, Gdir, 'montage');
title('Gradient Magnitude, Gmag (left), and Gradient Direction, Gdir (right), using Prewitt method')

%% Edge detection
I1_edge = edge(I1,'Canny');
figure; imshow(I1_edge);

%% Find pixel coordinates along a line

% point = round(siz/2);
point = [114 411];
angle = 20*pi/180;
range = 20;

if (angle>=-pi/4 && angle<=pi/4) || (angle<=-3*pi/4 && angle>=-pi) || ...
   (angle>=3*pi/4 && angle<=pi)
    % x is column and y is row
    if (angle>=-pi/4 && angle<0) || (angle>=3*pi/4 && angle<pi)
        xt = point(2):-1:point(2)-range;
        yt = round(tan(-angle)*(xt-point(2))+point(1));
        xb = point(2): 1:point(2)+range;
        yb = round(tan(-angle)*(xb-point(2))+point(1));   
    else
        xb = point(2):-1:point(2)-range;
        yb = round(tan(-angle)*(xb-point(2))+point(1));
        xt = point(2): 1:point(2)+range;
        yt = round(tan(-angle)*(xt-point(2))+point(1));   
    end
else
        yb = point(1): 1:point(1)+range;
        xb = round(((yb-point(1))/tan(-angle))+point(2));
        yt = point(1):-1:point(1)-range;
        xt = round(((yt-point(1))/tan(-angle))+point(2));
end

figure
hold on; imshow(I1_edge)
hold on;    plot(xb,yb,'r.')
hold on;    plot(xt,yt,'b.')


%% SSD Correlation along the line and evaluate optical flow

window = 15;        % should be odd number
best_match = [0 0];
best_corr  = Inf;
I1 = double(I1);
I2 = double(I2);

for i = 1:range+1
    
    % check if the search boundaries doesn't go out of the image
    if yb(i)-(window-1)/2>0 && yb(i)+(window-1)/2<=siz(1) && ...
            xb(i)-(window-1)/2>0 && xb(i)+(window-1)/2<=siz(2) && ...
            yt(i)-(window-1)/2>0 && yt(i)+(window-1)/2<=siz(1) && ...
            xt(i)-(window-1)/2>0 && xt(i)+(window-1)/2<=siz(2) && ...
            point(1)-(window-1)/2>0 && point(1)+(window-1)/2<=siz(2) && ...
            point(2)-(window-1)/2>0 && point(2)+(window-1)/2<=siz(2)
        
        corr_matr_l = I2(yb(i)-(window-1)/2:yb(i)+(window-1)/2,...
                         xb(i)-(window-1)/2:xb(i)+(window-1)/2)-...
                      I1(point(1)-(window-1)/2:point(1)+(window-1)/2,...
                         point(2)-(window-1)/2:point(2)+(window-1)/2);
        corr_l = sum(sum(corr_matr_l.^2));
        
        corr_matr_r = I2(yt(i)-(window-1)/2:yt(i)+(window-1)/2,...
                         xt(i)-(window-1)/2:xt(i)+(window-1)/2)-...
                      I1(point(1)-(window-1)/2:point(1)+(window-1)/2,...
                         point(2)-(window-1)/2:point(2)+(window-1)/2);
        corr_r = sum(sum(corr_matr_r.^2));
        
        % evaluate the lowest of correlation in both directions
        if corr_l < corr_r
            corr_d = corr_l;
            corr_match = [xb(i) yb(i)];
        else
            corr_d = corr_r;
            corr_match = [xt(i) yt(i)];
        end
        corr_d;
        % in case a better correlation is found
        if corr_d < best_corr
            best_corr  = corr_d
            best_match = corr_match;
        end
        
    end
end

best_corr
best_match
figure(1);     imshow(uint8(I2));
hold on;    plot(best_match(1),best_match(2),'bo');

figure(2);     imshow(uint8(I1));
hold on;    plot(point(2),point(1),'bo');
hold on;    plot(best_match(1),best_match(2),'bo');

flow_mag = dist(best_match,[point(2); point(1)])
flow_scale = 5;
hold on;
quiver(point(2),point(1),flow_scale*flow_mag*cos(angle),-flow_scale*...
    flow_mag*sin(angle),'Color','r');


%% Optical flow perpendicular to edges of image





