function [FOE, Flow_x, Flow_y] = optical_flow_edges(im1col, im2col, graphics)

%% Program to find the optical flow perpendicular to the edges in the image
%  Author: Ajith Anil Meera, 28th July 2017

if(~exist('graphics', 'var') || isempty(graphics))
    graphics = false;
end

%% Evaluate image gradient

% clc; close all; % better not to do this in a function
% im1col = imresize(im1col,.4);
% I2 = imresize(I2,0.4);

I1 = rgb2gray(im1col);
I2 = rgb2gray(im2col);

siz = size(I1);
% figure
% imshow(I)
[dx, dy] = imgradientxy(I1,'sobel');
Gdir = atan2(-dy,dx)*180/pi;    % in angles

% [Gmag, Gdir] = imgradient(I1,'prewitt');
% figure
% imshowpair(Gmag, Gdir, 'montage');
% title('Gradient Magnitude, Gmag (left), and Gradient Direction, Gdir (right), using Prewitt method')

%% Edge detection using Canny edge detector
I1_edge = edge(I1,'Canny',0.1);

% this can be sped up by using subsampling to get an idea of the mean abs
% gradients.
adx = abs(dx);
ady = abs(dy);
madx = mean(mean(adx));
mady = mean(mean(ady));
edge_factor = 2.5;
ADX = adx >= edge_factor * madx;
ADY = ady >= edge_factor * mady;
Sobel_edge = ADX | ADY;

if(graphics)
    figure();
    subplot(1,2,1);
    imshow(I1_edge);
    title('Canny edge');
    subplot(1,2,2);
    imshow(Sobel_edge);
    title('Sobel edge');
end

% figure; imshow(I1_edge);
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
            
            % SSD Correlation along the line and evaluate optical flow
            
            window = 15;        % should be odd number
            best_match = [0 0];
            best_corr  = Inf;
            best_dir   = 1;     % positive direction
            all_corr = zeros(1,2*range+1);
            
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
                        best_corr  = corr_d;
                        best_match = corr_match;
                    end
                    
                end
            end
            % Plot the correlation 
%             tot_points
%             h1 = figure(1); plot(-range:range,all_corr,'r-'); 
%             xlabel('Search range perpendicular to edge pixel');
%             ylabel('SSD Correlation');
            
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
            if flow_mag(tot_points)>flow_thresh
                tot_points = tot_points - 1;
            end
            
        end
    end
end
flow_mag = flow_mag(1:tot_points);
angle = angle(1:tot_points);
point = point(1:tot_points,:);

%% Plot optical flow with a stride

stride = 1;
flow_scale = 1;
if graphics
    figure; imshow(im1col);
    hold on;
    quiver(point(1:stride:tot_points,2)',point(1:stride:tot_points,1)',...
        flow_mag(1:stride:tot_points).*cos(angle(1:stride:tot_points)),...
        -flow_mag(1:stride:tot_points).*sin(angle(1:stride:tot_points)),flow_scale,'Color','r');
end

%% Estimate the resultant flow direction

Flow_x = sum(flow_mag(find(flow_mag)).*cos(angle(find(flow_mag))));
Flow_y = sum(flow_mag(find(flow_mag)).*sin(angle(find(flow_mag))));
Flow_dir = atan2(Flow_y,Flow_x);
Flow_mag = hypot(Flow_x,Flow_y)/size(find(flow_mag),2) % average of all flow vectors
% Flow_mag = mode(flow_mag(find(flow_mag)))              % mode 

Flow_dir*180/pi

%%
% is_flow = find((round(flow_mag)==mode(round(flow_mag(find(flow_mag))))|round(flow_mag)==10)...
%                 &round(angle,1)~=0);
% is_flow = find(round(flow_mag)==11&round(angle,1)~=0);            
% is_flow = (round(flow_mag)<3&round(angle,1)~=0); 

% is_flow = find(round(angle*180/pi)>0&round(angle*180/pi)<90);

% plot((point(is_flow,2)-tan(angle(is_flow)'+pi/2).*point(is_flow,1))+...
%     tan(angle(is_flow)'+pi/2).*(1:510),1:510);
% 
% 
% FOE = pinv([-tan(angle(is_flow)+pi/2)' ones(size(angle(is_flow)))'])*...
%                 (point(is_flow,2)-tan(angle(is_flow)'+pi/2).*point(is_flow,1))

flow_mag_thresh = 2;

% keep low flow magnitude edge points and remove singularity points
is_flow = (round(flow_mag)<flow_mag_thresh&round(angle-pi/2,1)~=0&round(angle+pi/2,1)~=0);

% intersection of tangents to the edges at points with least flow magnitude
FOE = pinv([-tan(angle(is_flow))' ones(size(angle(is_flow)))'])*...
    (point(is_flow,2)-tan(angle(is_flow)').*point(is_flow,1));

if graphics
    hold on;
    plot((point(is_flow,2)-tan(angle(is_flow)').*point(is_flow,1))+...
        tan(angle(is_flow)').*(1:510),1:510);
    hold on; plot(FOE(2),FOE(1),'g*');
    hold on; quiver(FOE(2),FOE(1),Flow_mag*cos(Flow_dir),-Flow_mag*sin(Flow_dir),10,'Color','k');
end


%% Evaluate if it a translation or scaling 

if(Flow_mag<2)
    disp('Pure scaling');
else 
    disp('Translation or scaling');
end



end
