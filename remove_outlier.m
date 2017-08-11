function best_inliers = remove_outlier(point, angle, is_flow)
%% Function to remove the outlier tangent lines at the edges with zero flow magnitude
%  Author: Ajith Anil Meera, 11th Aug 2017

%% RANSAC to remove outlier lines

samp_perc = .5;
n_inliers = round(size(is_flow(is_flow==1),2)*samp_perc);
if n_inliers<4
    n_inliers = size(find(is_flow),2);
end
best_inliers = zeros(1,n_inliers);
best_FOE = [0 0];
best_inters = Inf;

for kk = 1:10000
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

% remove lines with outlier angles from the result
std_thresh = 2; % sigma
best_inliers = best_inliers(angle(best_inliers)>=mean(angle(is_flow))-...
    std_thresh*std(angle(is_flow))&angle(best_inliers)<=mean(angle(is_flow))...
    +std_thresh*std(angle(is_flow)));

end