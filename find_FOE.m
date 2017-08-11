function FOE = find_FOE(point, flow_mag, angle,graphics)
%% Function to find the FOE, given the edge points, flow magnitude and direction
%  Author: Ajith Anil Meera, 11th Aug 2017

if_outlier = 1;
global frame;
filename = 'testAnimated.gif';

%% Estimate the resultant flow direction

Flow_x = sum(flow_mag(find(flow_mag)).*cos(angle(find(flow_mag))));
Flow_y = sum(flow_mag(find(flow_mag)).*sin(angle(find(flow_mag))));
Flow_dir = atan2(Flow_y,Flow_x);
Flow_mag = hypot(Flow_x,Flow_y)/size(find(flow_mag),2); % average of all flow vectors
% Flow_mag = mode(flow_mag(find(flow_mag)))              % mode 

% fprintf('Flow magnitude: %f',Flow_mag);
% fprintf('\nFlow direction(deg): %f',Flow_dir*180/pi);

flow_mag_thresh = 1;

% keep low flow magnitude edge points and remove singularity points
is_flow = (round(flow_mag)<flow_mag_thresh&round(angle-pi/2,1)~=0&round(angle+pi/2,1)~=0);

if if_outlier
    is_flow = remove_outlier(point,angle,is_flow);
end

% intersection of tangents to the edges at points with least flow magnitude
FOE = pinv([-tan(angle(is_flow))' ones(size(angle(is_flow)))'])*...
    (point(is_flow,2)-tan(angle(is_flow)').*point(is_flow,1));

if graphics
    figure(3); hold on;
%     plot((point(is_flow,2)-tan(angle(is_flow)').*point(is_flow,1))+...
%         tan(angle(is_flow)').*(1:510),1:510);
    figure(3); hold on; plot(FOE(2),FOE(1),'g*');
    figure(3); hold on; quiver(FOE(2),FOE(1),3*Flow_mag*cos(Flow_dir),...
        -3*Flow_mag*sin(Flow_dir),10,'Color','b');
end

% Evaluate if it a translation or scaling 
if std(angle(is_flow))*180/pi<10           % if most of the angles are almost parallel
    figure(3); hold on; text(50,50,'Pure translation','Color','g')
elseif(Flow_mag<1)                          % if flow is not sufficient to set all pixels to motion
    figure(3); hold on; text(50,50,'Pure scaling','Color','g');
else 
    figure(3); hold on; text(50,50,'Translation and scaling','Color','g');
end
h1=figure(3);
% Capture the figure and save as gif
frame_gif = getframe(h1);
im_gif = frame2im(frame_gif);
[imind_gif,cm_gif] = rgb2ind(im_gif,256);

% Write to the GIF File
if frame == 75
    imwrite(imind_gif,cm_gif,filename,'gif', 'Loopcount',inf);
else
    imwrite(imind_gif,cm_gif,filename,'gif','WriteMode','append');
end


end