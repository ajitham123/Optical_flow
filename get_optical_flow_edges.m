function [point, flow_x, flow_y] = get_optical_flow_edges(I1, I2, graphics)

if(~exist('graphics', 'var') || isempty(graphics))
    graphics = false;
end

%% Program to find the optical flow perpendicular to the edges in the image
%  Author: Ajith Anil Meera, 28th July 2017

%% Evaluate image gradient
siz = size(I1);
% figure
% imshow(I)
[Gmag, Gdir] = imgradient(I1,'prewitt');

if(graphics)
    figure
    imshowpair(Gmag, Gdir, 'montage');
    title('Gradient Magnitude, Gmag (left), and Gradient Direction, Gdir (right), using Prewitt method')
end

%% Edge detection using Canny edge detector
I1_edge = edge(I1,'Canny');
if(graphics)
    figure; imshow(I1_edge);
end
I1 = double(I1);
I2 = double(I2);

%% Evaluate Optical flow perpendicular to all edges in the image

range = 20;
tot_points = 0;
flow_mag = zeros(1,prod(siz));
angle = zeros(1,prod(siz));
point = zeros(prod(siz),2);
flow_thresh = 100;

for j=1:siz(1)
    for k=1:siz(2)
        if I1_edge(j,k)==1      % for all edges
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
                    xt = point(tot_points,2)-1:-1:point(tot_points,2)-range;
                    yt = round(tan(-angle(tot_points))*(xt-point(tot_points,2))+point(tot_points,1));
                    xb = point(tot_points,2)+1: 1:point(tot_points,2)+range;
                    yb = round(tan(-angle(tot_points))*(xb-point(tot_points,2))+point(tot_points,1));
                else
                    xb = point(tot_points,2)-1:-1:point(tot_points,2)-range;
                    yb = round(tan(-angle(tot_points))*(xb-point(tot_points,2))+point(tot_points,1));
                    xt = point(tot_points,2)+1: 1:point(tot_points,2)+range;
                    yt = round(tan(-angle(tot_points))*(xt-point(tot_points,2))+point(tot_points,1));
                end
            else
                yb = point(tot_points,1)+1: 1:point(tot_points,1)+range;
                xb = round(((yb-point(tot_points,1))/tan(-angle(tot_points)))+point(tot_points,2));
                yt = point(tot_points,1)-1:-1:point(tot_points,1)-range;
                xt = round(((yt-point(tot_points,1))/tan(-angle(tot_points)))+point(tot_points,2));
            end
            
            % SSD Correlation along the line and evaluate optical flow
            
            window = 15;        % should be odd number
            best_match = [0 0];
            best_corr  = Inf;
            best_dir   = 1;     % positive direction
            
            for i = 1:range
                
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
            if flow_mag(tot_points)>flow_thresh || isinf(best_corr) 
                tot_points = tot_points - 1;
            end
            
        end
    end
end
flow_mag = flow_mag(1:tot_points);
angle = angle(1:tot_points);
point = point(1:tot_points,:);
flow_x = flow_mag(1:tot_points).*cos(angle(1:tot_points));
flow_y = -flow_mag(1:tot_points).*sin(angle(1:tot_points));

%% Plot optical flow with a stride

if(graphics)
    stride = 2;
    flow_scale = .5;
    figure; imshow(I1 ./ 255);
    hold on;
    quiver(point(1:stride:tot_points,2)',point(1:stride:tot_points,1)',...
        flow_mag(1:stride:tot_points).*cos(angle(1:stride:tot_points)),...
        -flow_mag(1:stride:tot_points).*sin(angle(1:stride:tot_points)),flow_scale,'Color','r');
end

