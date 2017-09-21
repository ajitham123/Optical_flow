function best_inliers = remove_outlier(point, angle, is_flow)
%% Function to remove the outlier tangent lines at the edges with zero flow magnitude
%  Author: Ajith Anil Meera, 11th Aug 2017

%% RANSAC to remove outlier lines

samp_perc = 1;
n_inliers = round(size(is_flow(is_flow==1),2)*samp_perc);
if n_inliers<5
    disp('Too few inlier lines');
    n_inliers = size(find(is_flow),2);
end
% n_inliers = 30;
best_inliers = zeros(1,n_inliers);
best_inliers = find(is_flow);
best_FOE = [0 0];
best_inters = Inf;


% for kk = 1:1 % 20000
%     inliers = randsample(find(is_flow),n_inliers);
%     FOE = pinv([-tan(angle(inliers))' ones(size(angle(inliers)))'])*...
%         (point(inliers,2)-tan(angle(inliers)').*point(inliers,1));
%     error_inters = sum(abs(tan(angle(inliers))'*FOE(1)-FOE(2)+point(inliers,2)-...
%         tan(angle(inliers))'.*point(inliers,1)));
%     if error_inters < best_inters
%        best_inliers = inliers;
%        best_FOE = FOE;
%        best_inters = error_inters;
%     end
% end
% error_inters;
% size(best_inliers,2);
% best_FOE;

i = 1;
err_std_thresh = 2.4; % 2.4
error_distances = abs(tan(angle(best_inliers))'*best_FOE(1)-best_FOE(2)+...
        point(best_inliers,2)-tan(angle(best_inliers))'.*point(best_inliers,1))./sqrt(1+(tan(angle(best_inliers))').^2);
while mean(error_distances)>3 && i<70 % 3 and 70
% for i=1:10
    best_inliers = best_inliers(error_distances >= mean(error_distances)-err_std_thresh*...
        std(error_distances) & error_distances <= mean(error_distances)+err_std_thresh*...
        std(error_distances));
    best_FOE = pinv([-tan(angle(best_inliers))' ones(size(angle(best_inliers)))'])*...
        (point(best_inliers,2)-tan(angle(best_inliers)').*point(best_inliers,1));
    error_distances = abs(tan(angle(best_inliers))'*best_FOE(1)-best_FOE(2)+...
        point(best_inliers,2)-tan(angle(best_inliers))'.*point(best_inliers,1))./sqrt(1+(tan(angle(best_inliers))').^2);
    i = i+1;
%     err_std_thresh = err_std_thresh + i*0.1;
    
%     figure(3); imshow(imread('images/indoor.jpg')); hold on;
%     plot((point(best_inliers,2)-tan(angle(best_inliers)').*point(best_inliers,1))+...
%         tan(angle(best_inliers)').*(1:510),1:510);
%     figure(3); hold on; plot(point(best_inliers,2),point(best_inliers,1),'c*'); 
%     figure(3); hold on; plot(best_FOE(2),best_FOE(1),'g*');
%     pause(.002);
    
end
mean(error_distances);
size(best_inliers,2);

bbest_inters = Inf;
for kk = 1:500 % 20000
    inliers = randsample(best_inliers,200);
    FOE = pinv([-tan(angle(inliers))' ones(size(angle(inliers)))'])*...
        (point(inliers,2)-tan(angle(inliers)').*point(inliers,1));
    error_inters = sum(abs(tan(angle(inliers))'*FOE(1)-FOE(2)+point(inliers,2)-...
        tan(angle(inliers))'.*point(inliers,1))./sqrt(1+(tan(angle(inliers))').^2));
    if error_inters < bbest_inters
       bbest_inliers = inliers;
       bbest_inters = error_inters;
    end
end
best_inliers = bbest_inliers;


% remove lines with outlier angles from the result
% std_thresh = 1.5; % sigma
% best_inliers = best_inliers(angle(best_inliers)>=mean(angle(is_flow))-...
%     std_thresh*std(angle(is_flow))&angle(best_inliers)<=mean(angle(is_flow))...
%     +std_thresh*std(angle(is_flow)));

end