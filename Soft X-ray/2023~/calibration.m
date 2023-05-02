clear
date = 230118;
CalibrationPath = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/SXR_Images/',num2str(date),'/PositionCheck.tif');
CalibrationImage = imread(CalibrationPath);
[centers,radii]=FindFibers(CalibrationImage);
IW = round(min(radii));
centers = round(centers);

Center = zeros(2,8,2);
TimeRapsImage = zeros(2,8,2*IW,2*IW);

Center(1,:,:) = centers([1,3,6,8,9,11,14,16],:);
Center(2,:,:) = centers([2,4,5,7,10,12,13,15],:);

RawImage = CalibrationImage;
% disp(Center);
% whos RawImage
% whos TimeRapsImage

for i = 1:8
    UpRangeV = Center(1,i,1)-IW+1:Center(1,i,1)+IW;
    UpRangeH = Center(1,i,2)-IW+1:Center(1,i,2)+IW;
    DownRangeV = Center(2,i,1)-IW+1:Center(2,i,1)+IW;
    DownRangeH = Center(2,i,2)-IW+1:Center(2,i,2)+IW;
    TimeRapsImage(1,i,:,:) = RawImage(UpRangeV,UpRangeH,1);
    TimeRapsImage(2,i,:,:) = RawImage(DownRangeV,DownRangeH,1);
    TimeRapsImage(1,i,:,:) = rot90(TimeRapsImage(1,i,:,:),2);
    TimeRapsImage(2,i,:,:) = rot90(TimeRapsImage(2,i,:,:),2);
end

BackGround = cast(RawImage(1050-IW+1:1050+IW,120-IW+1:120+IW,1),'double');
BackGround = ones(size(BackGround))*mean(BackGround,'all');

N_projection = 80;
k = FindCircle(N_projection/2);
VectorImages = zeros(2,8,numel(k));
% vector_images_filtered = zeros(2,8,numel(k));

check_flag = true;

if check_flag
    f1=figure;
%     f.WindowState = 'maximized';
    f1.Position = [900,200,500,500];
end

resolution = N_projection/(IW*2);
measurement_error = zeros(2,8);

for i=1:8
    SingleImage1 = squeeze(TimeRapsImage(1,i,:,:))-BackGround;
    SingleImage2 = squeeze(TimeRapsImage(2,i,:,:))-BackGround;
    GrayImage1 = fliplr(SingleImage1(:,:));
    GrayImage2 = fliplr(SingleImage2(:,:));
    GrayImage1(GrayImage1<0) = 0;
    GrayImage2(GrayImage2<0) = 0;
    RoughImage1 = imresize(GrayImage1, resolution, 'nearest');
    RoughImage1 = cast(RoughImage1,'double');
    RoughImage2 = imresize(GrayImage2, resolution, 'nearest');
    RoughImage2 = cast(RoughImage2,'double');
    
%     whos RoughImage1
%     whos SD
%     
%     RoughImage1 = RoughImage1./fliplr(SD./100);
%     RoughImage2 = RoughImage2./fliplr(SD./100);
    
    if check_flag
        figure(f1);
        i_str = num2str(i);
        title1 = strcat('1,',i_str);
        title2 = strcat('2,',i_str);
        subplot(4,4,2*(i-1)+1);imagesc(RoughImage1);title(title1);
        caxis([50,60]);
        subplot(4,4,2*(i-1)+2);imagesc(RoughImage2);title(title2);
        caxis([50,60]);
%         if i == 8
%             colorbar;
%         end
    end

    VectorImages(1,i,:) = RoughImage1(k);
    VectorImages(2,i,:) = RoughImage2(k);
    measurement_error(1,i) = std(VectorImages(1,i,:));
    measurement_error(2,i) = std(VectorImages(2,i,:));
%     rough_images_filtered1 = imgaussfilt(RoughImage1);
%     rough_images_filtered2 = imgaussfilt(RoughImage2);
%     vector_images_filtered(1,i,:) = rough_images_filtered1(k);
%     vector_images_filtered(2,i,:) = rough_images_filtered2(k);
end

% mean_filtered = mean(vector_images_filtered,'all');
% calibration_factor_filtered = mean_filtered./vector_images_filtered;
MeanIntensity = mean(VectorImages,'all');
CalibrationFactor = MeanIntensity./VectorImages;

% calibrated_images = calibration_factor_filtered.*VectorImages;

% if check_flag
%     f2=figure;
% %     f.WindowState = 'maximized';
%     f2.Position = [900,200,500,500];
% end
% 
% for i=1:8
%     SingleImage1 = squeeze(TimeRapsImage(1,i,:,:))-BackGround;
%     SingleImage2 = squeeze(TimeRapsImage(2,i,:,:))-BackGround;
%     GrayImage1 = fliplr(SingleImage1(:,:));
%     GrayImage2 = fliplr(SingleImage2(:,:));
%     GrayImage1(GrayImage1<0) = 0;
%     GrayImage2(GrayImage2<0) = 0;
%     RoughImage1 = imresize(GrayImage1, resolution, 'nearest');
%     RoughImage1 = cast(RoughImage1,'double');
%     RoughImage2 = imresize(GrayImage2, resolution, 'nearest');
%     RoughImage2 = cast(RoughImage2,'double');
%     
% %     whos RoughImage1
% %     whos CalibrationFactor
% %     A = RoughImage1(k);
% %     B = CalibrationFactor(1,i,:);
% %     whos A
% %     whos B
% %     A = RoughImage1(k).*CalibrationFactor(1,i,:);
% %     whos A
%     RoughImage1(k) = RoughImage1(k).*squeeze(CalibrationFactor(1,i,:));
%     RoughImage2(k) = RoughImage2(k).*squeeze(CalibrationFactor(2,i,:));
%     
% %     whos RoughImage1
% %     whos SD
% %     
% %     RoughImage1 = RoughImage1./fliplr(SD./100);
% %     RoughImage2 = RoughImage2./fliplr(SD./100);
%     
%     if check_flag
%         figure(f2);
%         i_str = num2str(i);
%         title1 = strcat('1,',i_str);
%         title2 = strcat('2,',i_str);
%         subplot(4,4,2*(i-1)+1);imagesc(RoughImage1);title(title1);
%         caxis([50,60]);
%         subplot(4,4,2*(i-1)+2);imagesc(RoughImage2);title(title2);
%         caxis([50,60]);
% %         if i == 8
% %             colorbar;
% %         end
%     end
% end

function k = FindCircle(L)
R = zeros(2*L);
for i = 1:2*L
    for j = 1:2*L
        R(i,j) = sqrt((L-i+0.5)^2+(j-L-0.5)^2);
    end
end
% figure;imagesc(R)
k = find(R<L);
end



