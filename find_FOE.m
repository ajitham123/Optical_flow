function FOE = find_FOE(point, flow_mag, angle,graphics)
%% Function to find the FOE, given the edge points, flow magnitude and direction
%  Author: Ajith Anil Meera, 11th Aug 2017

if_outlier = 1;
%% Estimate the resultant flow direction

Flow_x = sum(flow_mag(find(flow_mag)).*cos(angle(find(flow_mag))));
Flow_y = sum(flow_mag(find(flow_mag)).*sin(angle(find(flow_mag))));
Flow_dir = atan2(Flow_y,Flow_x);
Flow_mag = hypot(Flow_x,Flow_y)/size(find(flow_mag),2); % average of all flow vectors
% Flow_mag = mode(flow_mag(find(flow_mag)))              % mode 

fprintf('Flow magnitude: %f',Flow_mag);
fprintf('\nFlow direction(deg): %f',Flow_dir*180/pi);

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
    hold on;
    plot((point(is_flow,2)-tan(angle(is_flow)').*point(is_flow,1))+...
        tan(angle(is_flow)').*(1:510),1:510);
    hold on; plot(FOE(2),FOE(1),'g*');
    hold on; quiver(FOE(2),FOE(1),Flow_mag*cos(Flow_dir),-Flow_mag*sin(Flow_dir),10,'Color','k');
end

% Evaluate if it a translation or scaling 
if std(angle(is_flow))*180/pi<7.5
    fprintf('\nPure translation')
elseif(Flow_mag<2)
    fprintf('\nPure scaling');
else 
    fprintf('\nTranslation or scaling or both');
end

end