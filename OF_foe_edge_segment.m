%% Program to find the optical flow perpendicular to the edges of all frames
%  in a video sequence
%  Author: Ajith Anil Meera, 2nd August 2017

%% Evaluate image gradient

clear all; clc; close all;
im_scale = 0.4;     % scale down for high resolution images
window = 15;        % should be odd number, window for SSD correlation

filename = 'testAnimated.gif';
for fram=25:75
    im1col = imread(['C:\Users\Ajith A M\OneDrive\Windows desktop backup\MAVlab\frames\frames',int2str(fram),'.jpg']);
    I2 = rgb2gray(imread(['C:\Users\Ajith A M\OneDrive\Windows desktop backup\MAVlab\frames\frames',int2str(fram+1),'.jpg']));
    im1col = imresize(im1col,im_scale);
    I2 = imresize(I2,im_scale);
    I1 = rgb2gray(im1col);
    
    siz = size(I1);
    [dx, dy] = imgradientxy(I1,'sobel');
    Gdir = atan2(-dy,dx)*180/pi;    % in angles
    %[~, Gdir] = imgradient(I1,'prewitt');

    % Edge detection using Canny edge detector
    I1_edge = edge(I1,'Canny',.3);
    % figure; imshow(I1_edge);
    I1 = double(I1);
    I2 = double(I2);
    
    % Evaluate Optical flow perpendicular to all edges in the image
    
    range = 20;
    tot_points = 0;
    prev_tot_points = 0;
    n_edges = 3;
    flow_mag = zeros(1,prod(siz));
    angle = zeros(1,prod(siz));
    point = zeros(prod(siz),2);
    flow_thresh = sqrt(2)*range;
    FOE = zeros(n_edges,2);
    
    % Dilate and connect weak edges. Label the connected edges
    I1_edge_label = bwlabel(imdilate(I1_edge,strel('line',3,3)),4);
    
    % Evaluate largest connected edges
    edge_count = hist(I1_edge_label(I1_edge_label>0),unique(I1_edge_label(I1_edge_label>0)));
    [~, edge_index] = sort(edge_count,'descend');
    
    for edd = 1:n_edges
        
        % Segment the edge out for processing based on the edge index
        I1_edge_seg = (I1_edge_label==edge_index(edd));
        
        for j=ceil(window/2):siz(1)-floor(window/2)
            for k=ceil(window/2):siz(2)-floor(window/2)
                % Boundary of image where correlation cannot be calculated due to the
                % window size is left out from flow estimation
                
                if I1_edge_seg(j,k)==1      % for all points in the segmented edge
                    tot_points = tot_points+1;
                    angle(tot_points) = Gdir(j,k)*pi/180;
                    point(tot_points,:) = [j k];      % image coordinates
                    
                    % find search points along the gradient direction - all points
                    % found will fall inside a cube of size 'range' centred about
                    % 'point'
                    if (angle(tot_points)>=-pi/4 && angle(tot_points)<=pi/4)  || ...
                            (angle(tot_points)<=-3*pi/4 && angle(tot_points)>=-pi) || ...
                            (angle(tot_points)>=3*pi/4 && angle(tot_points)<=pi)
                        
                        if (angle(tot_points)>=-pi/4 && angle(tot_points)<0) || ...
                                (angle(tot_points)>=3*pi/4 && angle(tot_points)<pi)
                            % x is column and y is row
                            xt = point(tot_points,2):-1:point(tot_points,2)-range;
                            yt = round(tan(-angle(tot_points))*(xt-point(tot_points,2))+point(tot_points,1));
                            xb = point(tot_points,2): 1:point(tot_points,2)+range;
                            yb = round(tan(-angle(tot_points))*(xb-point(tot_points,2))+point(tot_points,1));
                        else
                            xb = point(tot_points,2):-1:point(tot_points,2)-range;
                            yb = round(tan(-angle(tot_points))*(xb-point(tot_points,2))+point(tot_points,1));
                            xt = point(tot_points,2): 1:point(tot_points,2)+range;
                            yt = round(tan(-angle(tot_points))*(xt-point(tot_points,2))+point(tot_points,1));
                        end
                    else
                        yb = point(tot_points,1): 1:point(tot_points,1)+range;
                        xb = round(((yb-point(tot_points,1))/tan(-angle(tot_points)))+point(tot_points,2));
                        yt = point(tot_points,1):-1:point(tot_points,1)-range;
                        xt = round(((yt-point(tot_points,1))/tan(-angle(tot_points)))+point(tot_points,2));
                    end
                    
                    % SSD Correlation along the line and evaluate optical flow
                    best_match = [0 0];
                    best_corr  = Inf;
                    best_dir   = 1;     % positive direction
                    
                    for i = 1:range+1
                        
                        % check if the search boundaries doesn't go out of the image
                        if yb(i)-(window-1)/2>0 && yb(i)+(window-1)/2<=siz(1) && ...
                                xb(i)-(window-1)/2>0 && xb(i)+(window-1)/2<=siz(2) && ...
                                yt(i)-(window-1)/2>0 && yt(i)+(window-1)/2<=siz(1) && ...
                                xt(i)-(window-1)/2>0 && xt(i)+(window-1)/2<=siz(2) && ...
                                point(tot_points,1)-(window-1)/2>0 && point(tot_points,1)+(window-1)/2<=siz(2) && ...
                                point(tot_points,2)-(window-1)/2>0 && point(tot_points,2)+(window-1)/2<=siz(2)
                            
                            corr_matr_b = I2(yb(i)-(window-1)/2:yb(i)+(window-1)/2,...
                                xb(i)-(window-1)/2:xb(i)+(window-1)/2)-...
                                I1(point(tot_points,1)-(window-1)/2:point(tot_points,1)+(window-1)/2,...
                                point(tot_points,2)-(window-1)/2:point(tot_points,2)+(window-1)/2);
                            corr_b = sum(sum(corr_matr_b.^2));
                            
                            corr_matr_t = I2(yt(i)-(window-1)/2:yt(i)+(window-1)/2,...
                                xt(i)-(window-1)/2:xt(i)+(window-1)/2)-...
                                I1(point(tot_points,1)-(window-1)/2:point(tot_points,1)+(window-1)/2,...
                                point(tot_points,2)-(window-1)/2:point(tot_points,2)+(window-1)/2);
                            corr_t = sum(sum(corr_matr_t.^2));
                            
                            % evaluate the lowest of correlation in both directions
                            % along the gradient
                            if corr_b < corr_t
                                corr_dir = 0;       % negative direction
                                corr_d = corr_b;
                                corr_match = [xb(i) yb(i)];
                            else
                                corr_dir = 1;       % flag for top portion (positive direction)
                                corr_d = corr_t;    % correlation value in the direction
                                corr_match = [xt(i) yt(i)];
                            end
                            
                            % in case a better correlation is found update best
                            % matches
                            if corr_d < best_corr
                                best_dir   = corr_dir;
                                best_corr  = corr_d;
                                best_match = corr_match;
                            end
                            
                        end
                    end
                    
                    % If best correlation is in the opposite direction of gradient,
                    % update the angle
                    if best_dir == 1 && angle(tot_points)<0
                        angle(tot_points) = angle(tot_points)+pi;
                    elseif best_dir == 0 && angle(tot_points)>=0
                        angle(tot_points) = angle(tot_points)-pi;
                    end
                    
                    flow_mag(tot_points) = dist(best_match,[point(tot_points,2); point(tot_points,1)]);
                    
                    % in case flow magnitude is too large - particularly near image
                    % boundary - remove the flow
                    if flow_mag(tot_points)>flow_thresh
                        tot_points = tot_points - 1;
                    end
                    
                end
            end
        end
        
        
        % Estimate Focus of expansion by the assumption that the centre of all
        % non zero flow points found at the edges in the frame is the required
        % FOE. The points are assumed to be in a circle.
        
        [FOE_x,FOE_y] = circfit(point(prev_tot_points+find(flow_mag(prev_tot_points+1:tot_points)),1)...
                               ,point(prev_tot_points+find(flow_mag(prev_tot_points+1:tot_points)),2));
        FOE(edd,:) = [FOE_x FOE_y];
        
        % FoE by minimising least square error for the set of optical flow
        % vectors
%         is_flow = prev_tot_points+find(flow_mag(prev_tot_points+1:tot_points));
%         FOE(edd,:) = pinv([-tan(angle(is_flow)+pi/2)' ones(size(angle(is_flow)))'])*...
%             (point(is_flow,2)-tan(angle(is_flow)'+pi/2).*point(is_flow,1));
        
        prev_tot_points = tot_points;
    end
    flow_mag = flow_mag(1:tot_points);
    angle = angle(1:tot_points);
    point = point(1:tot_points,:);
    
    % Plot optical flow with a stride
    
    stride = 2;
    flow_scale = 1;
    h1 = figure(1); imshow(im1col);
    hold on;
    quiver(point(1:stride:tot_points,2)',point(1:stride:tot_points,1)',...
        flow_mag(1:stride:tot_points).*cos(angle(1:stride:tot_points)),...
        -flow_mag(1:stride:tot_points).*sin(angle(1:stride:tot_points)),flow_scale,'Color','b');
%     hold on;
%     plot((point(is_flow,2)-tan(angle(is_flow)'+pi/2).*point(is_flow,1))+...
%         tan(angle(is_flow)'+pi/2).*(1:510),1:510);
    
    hold on;
    plot(FOE(:,2), FOE(:,1),'g*'); hold on;
    pause(.0001)
    
    % Capture the figure and save as gif
    frame_gif = getframe(h1);
    im_gif = frame2im(frame_gif);
    [imind_gif,cm_gif] = rgb2ind(im_gif,256);
    
    % Write to the GIF File
    if fram == 1
        imwrite(imind_gif,cm_gif,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind_gif,cm_gif,filename,'gif','WriteMode','append');
    end
end





