function [point, flow_mag, angle] = get_optical_flow_edges(I1, I2, graphics,if_sub_pixel,if_ssd)
%% Program to find the optical flow perpendicular to the edges in the image
%  Author: Ajith Anil Meera, 28th July 2017
% figure(3); imshow(I1);
%     hold on;
if(~exist('graphics', 'var') || isempty(graphics))
    graphics = false;
    if_sub_pixel = false;
    if_ssd = true;
end

%% Evaluate image gradient
if(size(I1,3) == 3)
    I1 = rgb2gray(I1);
    I2 = rgb2gray(I2);
end

%  I1 = imgaussfilt(I1,2);
%  I2 = imgaussfilt(I2,2);
% I1 = imboxfilt(I1,3);
% I1 = medfilt2(I1);
siz = size(I1);
[dx, dy] = imgradientxy(I1,'sobel');
Gmag = sqrt(dx.^2 + dy.^2);
Gdir = atan2(-dy,dx)*180/pi;    % in angles
% dx(201,1:500)
% dy(201,1:500)
% Gdir(201,1:500)


% if(graphics)
%     figure
%     imshowpair(Gmag, Gdir, 'montage');
%     title('Gradient Magnitude, Gmag (left), and Gradient Direction, Gdir (right), using Sobel method')
% end

%% Edge detection using Sobel filter

% this can be sped up by using subsampling to get an idea of the mean abs
% gradients.
if_canny = 0;

if if_canny
    I1_edge = edge(I1,'Canny',0.1);
else
    adx = abs(dx);
    ady = abs(dy);
    madx = mean(mean(adx));
    mady = mean(mean(ady));
    edge_factor = 4;
    ADX = adx >= edge_factor * madx;
    ADY = ady >= edge_factor * mady;
    I1_edge = ADX | ADY;
%     I1_edge = imerode(I1_edge,strel('line',2,90));
end


% if(graphics)
%     figure(2); imshow(I1_edge);
% end
I1 = double(I1);
I2 = double(I2);

%% Evaluate Optical flow perpendicular to all edges in the image

range = 20;
tot_points = 0;
flow_mag = zeros(1,prod(siz));
angle = zeros(1,prod(siz));
point = zeros(prod(siz),2);
flow_thresh = sqrt(2)*range;
filename = 'testAnimated.gif';

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
%             xb
%             pause(1)

            % SSD Correlation along the line and evaluate optical flow
            
            window = 15;        % should be odd number
            best_match = [0 0];
            best_corr  = Inf;
            best_dir   = 1;     % positive direction
            best_i     = 0;     % best index
            all_corr = zeros(1,2*range+1);
            
            for i = 1:range+1
                
                % check if the search boundaries doesn't go out of the image
                if yb(i)-(window-1)/2>0 && yb(i)+(window-1)/2<=siz(1) && ...
                        xb(i)-(window-1)/2>0 && xb(i)+(window-1)/2<=siz(2) && ...
                        yt(i)-(window-1)/2>0 && yt(i)+(window-1)/2<=siz(1) && ...
                        xt(i)-(window-1)/2>0 && xt(i)+(window-1)/2<=siz(2) && ...
                        point(tot_points,1)-(window-1)/2>4 && point(tot_points,1)+(window-1)/2<=siz(1)-4 && ...
                        point(tot_points,2)-(window-1)/2>4 && point(tot_points,2)+(window-1)/2<=siz(2)-4 
                        % '4' above kept as tolerance against false zero flow
                        % detection near image boundary
                    
                        if if_ssd   
                            % SSD correlation
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
                            
                        else
                            % SAD corelation
                            corr_matr_b = I2(yb(i)-(window-1)/2:yb(i)+(window-1)/2,...
                                xb(i)-(window-1)/2:xb(i)+(window-1)/2)-...
                                I1(point(tot_points,1)-(window-1)/2:point(tot_points,1)+(window-1)/2,...
                                point(tot_points,2)-(window-1)/2:point(tot_points,2)+(window-1)/2);
                            corr_b = sum(sum(abs(corr_matr_b)));
                            
                            corr_matr_t = I2(yt(i)-(window-1)/2:yt(i)+(window-1)/2,...
                                xt(i)-(window-1)/2:xt(i)+(window-1)/2)-...
                                I1(point(tot_points,1)-(window-1)/2:point(tot_points,1)+(window-1)/2,...
                                point(tot_points,2)-(window-1)/2:point(tot_points,2)+(window-1)/2);
                            corr_t = sum(sum(abs(corr_matr_t)));
                        end
                    
                    % store all the correlation values to plot
                    all_corr(range+2-i) = corr_b;
                    all_corr(range+i) = corr_t;
                    
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
                        best_i     = i;
                        best_corr  = corr_d;
                        best_match = corr_match;
                    end
                    
                end
            end
            
            % sub pixel flow to refine flow magnitudes
            if if_sub_pixel
                if best_i-1>=1 && best_i+1<=range+1   % matrix dimension is correct
                    if best_dir == 0                  % if best match is at bottom portion
                        ii = range+2-best_i-1:range+2-best_i+1;
                        cc = all_corr(ii);
                        fit  = polyfit(ii,cc,2);      % fit a parabola
                        i_min = -fit(2)/(2*fit(1));   % minimum correlation index
                        i_incr = i_min-(range+2-best_i);    % gives the percent increment required
                        % perform an increment in best matches
                        if i_incr>0
                            x_min = best_match(1)+i_incr*(xb(best_i+1)-best_match(1));
                            y_min = best_match(2)+i_incr*(yb(best_i+1)-best_match(2));
                        else
                            x_min = best_match(1)+i_incr*(xb(best_i-1)-best_match(1));
                            y_min = best_match(2)+i_incr*(yb(best_i-1)-best_match(2));
                        end
                    else
                        ii = range+best_i-1:range+best_i+1;
                        cc = all_corr(ii);
                        fit  = polyfit(ii,cc,2);
                        i_min = -fit(2)/(2*fit(1));
                        i_incr = i_min-(range+best_i);
                        if i_incr>0
                            x_min = best_match(1)+i_incr*(xt(best_i+1)-best_match(1));
                            y_min = best_match(2)+i_incr*(yt(best_i+1)-best_match(2));
                        else
                            x_min = best_match(1)+i_incr*(xt(best_i-1)-best_match(1));
                            y_min = best_match(2)+i_incr*(yt(best_i-1)-best_match(2));
                        end
                    end
                    best_match = [x_min y_min];
                end
            end
            
            % Plot the correlation 
%             if tot_points == 8485    % 990, 8470
%                 point(tot_points,:)
%                 dist(best_match,[point(tot_points,2); point(tot_points,1)])
%                 h1 = figure(1); plot(-range:range,all_corr,'r-');
%                 xlabel('Search range perpendicular to edge pixel');
%                 ylabel('SSD Correlation'); pause(.00001);
%             end
            
            % Capture the figure and save as gif
%             frame_gif = getframe(h1);
%             im_gif = frame2im(frame_gif);
%             [imind_gif,cm_gif] = rgb2ind(im_gif,256);
%             
%             % Write to the GIF File
%             if tot_points == 200
%                 imwrite(imind_gif,cm_gif,filename,'gif', 'Loopcount',inf);
%             else
%                 imwrite(imind_gif,cm_gif,filename,'gif','WriteMode','append');
%             end
%             pause(.00001);

            
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

%% Plot optical flow with a stride

stride = 10;
flow_scale = 1;
if graphics
    figure(3); imshow(I1/255);
    hold on;
    quiver(point(1:stride:tot_points,2)',point(1:stride:tot_points,1)',...
        flow_mag(1:stride:tot_points).*cos(angle(1:stride:tot_points)),...
        -flow_mag(1:stride:tot_points).*sin(angle(1:stride:tot_points)),flow_scale,'Color','r');
end



end
