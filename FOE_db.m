function I2 = FOE_db(I1,FOE,scale)
%% Function to generate FOE image database, given scale and FOE
%  Author: Ajith Anil Meera, 26th August 2017

siz1 = size(I1);
% imshow(I1); hold on; plot(FOE(2),FOE(1),'r*');

% placing the FOE at image centre
siz2 = siz1(1:2)+2*abs(FOE-siz1(1:2)/2);
I2 = zeros([siz2 3]);
I2(siz2(1)/2-(FOE(1)-siz1(1)/2)+1-siz1(1)/2:siz2(1)/2-(FOE(1)-siz1(1)/2)+siz1(1)/2,...
    siz2(2)/2-(FOE(2)-siz1(2)/2)+1-siz1(2)/2:siz2(2)/2-(FOE(2)-siz1(2)/2)+siz1(2)/2,:)=I1 ;
I2 = uint8(I2);
% figure(2); imshow(I2); hold on; plot(siz2(2)/2,siz2(1)/2,'r*');

% scaling about FOE
I2 = imresize(I2,scale);
siz3 = size(I2);
centre = siz3(1:2)/2;
I2 = I2(floor(centre(1)-siz2(1)/2)+1:floor(centre(1)+siz2(1)/2),...
    floor(centre(2)-siz2(2)/2)+1:floor(centre(2)+siz2(2)/2),:);
% figure(3); imshow(I2); hold on;

% crop out to the initial image size
I2 = I2(siz2(1)/2-(FOE(1)-siz1(1)/2)+1-siz1(1)/2:siz2(1)/2-(FOE(1)-siz1(1)/2)+siz1(1)/2,...
    siz2(2)/2-(FOE(2)-siz1(2)/2)+1-siz1(2)/2:siz2(2)/2-(FOE(2)-siz1(2)/2)+siz1(2)/2,:);
% imshow(I2); hold on; plot(FOE(2),FOE(1),'r*');

% imwrite(I2,['images/indoor_foe_' int2str(FOE(1)) '_' int2str(FOE(2)) '.jpg']);

end
