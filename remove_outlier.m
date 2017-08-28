function best_inliers = remove_outlier(point, angle, is_flow)
%% Function to remove the outlier tangent lines at the edges with zero flow magnitude
%  Author: Ajith Anil Meera, 11th Aug 2017

%% RANSAC to remove outlier lines

% samp_perc = .05;
% n_inliers = round(size(is_flow(is_flow==1),2)*samp_perc);
% if n_inliers<5
%     disp('Too few inlier lines');
%     n_inliers = size(find(is_flow),2);
% end
n_inliers = 30;
best_inliers = zeros(1,n_inliers);
best_FOE = [0 0];
best_inters = Inf;

for kk = 1:20000
    inliers = randsample(find(is_flow),n_inliers);
    FOE = pinv([-tan(angle(inliers))' ones(size(angle(inliers)))'])*...
        (point(inliers,2)-tan(angle(inliers)').*point(inliers,1));
    error_inters = sum(abs(tan(angle(inliers))'*FOE(1)-FOE(2)+point(inliers,2)-...
        tan(angle(inliers))'.*point(inliers,1)));
    if error_inters < best_inters
       best_inliers = inliers;
       best_FOE = FOE;
       best_inters = error_inters;
    end
end
error_inters;
size(best_inliers,2);
best_FOE;

i = 1;
err_std_thresh = 1.5;
error_distances = abs(tan(angle(best_inliers))'*best_FOE(1)-best_FOE(2)+...
        point(best_inliers,2)-tan(angle(best_inliers))'.*point(best_inliers,1));
while mean(error_distances)>5 && i<2
% for i=1:10
    best_inliers = best_inliers(error_distances >= mean(error_distances)-err_std_thresh*...
        std(error_distances) & error_distances <= mean(error_distances)+err_std_thresh*...
        std(error_distances));
    best_FOE = pinv([-tan(angle(best_inliers))' ones(size(angle(best_inliers)))'])*...
        (point(best_inliers,2)-tan(angle(best_inliers)').*point(best_inliers,1));
    error_distances = abs(tan(angle(best_inliers))'*best_FOE(1)-best_FOE(2)+...
        point(best_inliers,2)-tan(angle(best_inliers))'.*point(best_inliers,1));
    i = i+1;
end
mean(error_distances);
size(best_inliers,2);


% remove lines with outlier angles from the result
std_thresh = 5; % sigma
best_inliers = best_inliers(angle(best_inliers)>=mean(angle(is_flow))-...
    std_thresh*std(angle(is_flow))&angle(best_inliers)<=mean(angle(is_flow))...
    +std_thresh*std(angle(is_flow)));

end