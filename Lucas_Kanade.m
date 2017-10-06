% clc; clear all; close all;
graphics = false;
stride = 50;
ww = 40;
w = round(ww/2);

% load images:
fr1 = imread('images/indoor.jpg');
h_map = zeros(floor(size(fr1,1)/stride),floor(size(fr1,2)/stride));

for ii=stride:stride:size(fr1,1)
    for jj=stride:stride:size(fr1,2)
        real_FOE = [ii,jj];
        fr2 = FOE_db(fr1,real_FOE,1.02);
        
        corn = detectHarrisFeatures(rgb2gray(fr1),'MinQuality',0.0000001,'ROI',...
            [w+1,w+1,size(fr1,2)-2*w-1,size(fr1,1)-2*w-1]);
        im1 = im2double(rgb2gray(fr1));
        im2 = im2double(rgb2gray(fr2));
        
        % Lucas Kanade Here
        % for each point, calculate I_x, I_y, I_t
        Ix_m = conv2(im1,[-1 1; -1 1], 'valid'); % partial on x
        Iy_m = conv2(im1, [-1 -1; 1 1], 'valid'); % partial on y
        It_m = conv2(im1, ones(2), 'valid') + conv2(im2, -ones(2), 'valid'); % partial on t
        u = zeros(size(im1));
        v = zeros(size(im2));
        point = zeros(corn.Count,2);
        angle = zeros(1,corn.Count);
        
        % within window ww * ww
        for k = 1:corn.Count    % w+1:size(Ix_m,1)-w
            %for j = w+1:size(Ix_m,2)-w
            
            i = round(corn.Location(k,2));
            j = round(corn.Location(k,1));
            point(k,:) = [i j];
            
            Ix = Ix_m(i-w:i+w, j-w:j+w);
            Iy = Iy_m(i-w:i+w, j-w:j+w);
            It = It_m(i-w:i+w, j-w:j+w);
            
            Ix = Ix(:);
            Iy = Iy(:);
            b = -It(:); % get b here
            
            A = [Ix Iy]; % get A here
            nu = pinv(A)*b; % get velocity here
            
            u(i,j)=nu(1);
            v(i,j)=nu(2);
            
            angle(k) = atan2(-v(i,j),u(i,j))+(pi/2);
        end
        
        % downsize u and v
        u_deci = u(1:1:end, 1:1:end);
        v_deci = v(1:1:end, 1:1:end);
        
        % get coordinate for u and v in the original frame
        [m, n] = size(im1);
        [X,Y] = meshgrid(1:n, 1:m);
        X_deci = X(1:1:end, 1:1:end);
        Y_deci = Y(1:1:end, 1:1:end);
        
        if graphics
            figure(3);    imshow(fr1);   hold on;
            % corn.plot; hold on
            plot(round(corn.Location(:,1)),round(corn.Location(:,2)),'b+');
            
            % draw the velocity vectors
            flow_scale = 100;
            quiver(X_deci, Y_deci, u_deci,v_deci, 'y')
        end
        
        % plot(point(300,2),point(300,1),'b*');
        % angle(300)*180/pi
        
        FOE = find_FOE(point,zeros(1,corn.Count),angle,graphics);
        h_map(ii/stride,jj/stride) = norm(real_FOE'-FOE);
        fprintf("\n i = %d, j = %d, FOE error = %f",ii,jj,h_map(ii/stride,jj/stride));
    end
end

fprintf("\n Mean error = %f ",mean(h_map(:)));
heatmap(h_map); 
figure; histogram(h_map,'BinWidth',10);
xlabel 'FOE error (in pixels)'
ylabel 'Count'

