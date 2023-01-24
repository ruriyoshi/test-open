function [VectorImage1,VectorImage2] = get_SXRImage(date,number,SXRfilename,filter)
% date = 230119;
% number = 5;
% SXRfilename = '/Users/shinjirotakeda/Library/CloudStorage/OneDrive-TheUniversityofTokyo/Documents/SXR_Images/230119/shot024.tif';
% filter = false;

% % 画像を切り取る
N_projection = 80;
% VectorImages = CutImage(date,shot,N_projection/230,false);
% CheckFlag = true;
CheckFlag = false;
VectorImages = CutImage(date,N_projection,CheckFlag,SXRfilename,filter);
VectorImages1 = squeeze(VectorImages(1,:,:));
VectorImages2 = squeeze(VectorImages(2,:,:));
VectorImage1 = VectorImages1(number,:);
VectorImage2 = VectorImages2(number,:);
end

function VectorImages = CutImage(date,N_projection,check_flag,SXRfilename,filter)

RawImage = imread(SXRfilename);

if filter
    figure;imagesc(RawImage);
    [RawImage,~] = imnlmfilt(RawImage,'SearchWindowSize',91,'ComparisonWindowSize',15);
    figure;imagesc(RawImage);
% else
%     RawImage = RawImage;
end

% IW = 115; %Image Width
% IW = 120; %Image Width

PositionPath = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/SXR_Images/',num2str(date),'/PositionCheck.tif');
CalibrationImage = imread(PositionPath);
[centers,radii]=FindFibers(CalibrationImage);
IW = round(mean(radii));
centers = round(centers);
IW = 80;

Center = zeros(2,8,2);
TimeRapsImage = zeros(2,8,2*IW,2*IW);

Center(1,:,:) = centers([1,3,6,8,9,11,14,16],:);
Center(2,:,:) = centers([2,4,5,7,10,12,13,15],:);

Center(2,1,:) = [530,1520];
Center(2,2,:) = [1060,1530];
Center(2,3,:) = [1070,1160];
Center(2,4,:) = [530,1150];
Center(2,5,:) = [530,770];
Center(2,6,:) = [1070,780];
Center(2,7,:) = [1070,410];
Center(2,8,:) = [530,380];
Center(1,:,1) = Center(2,:,1) - 190;
Center(1,:,2) = Center(2,:,2) + 5;

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

k = FindCircle(N_projection/2);
VectorImages = zeros(2,8,numel(k));

if check_flag
    f1=figure;
%     f.WindowState = 'maximized';
    f1.Position = [900,200,500,500];
    f2=figure;
%     f.WindowState = 'maximized';
    f2.Position = [900,200,500,500];
end

CalibrationPath = strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/SXR_Images/',num2str(date),'/CalibrationFactor.mat');
if exist(CalibrationPath,'file')
    load(CalibrationPath,'CalibrationFactor');
else
    CalibrationFactor = clc_CalibrationFactor(date,N_projection);
end

resolution = N_projection/(IW*2);

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
%         caxis([50,60]);
        subplot(4,4,2*(i-1)+2);imagesc(RoughImage2);title(title2);
%         caxis([50,60]);
%         if i == 8
%             colorbar;
%         end
        figure(f2);
        RoughCalibrated1 = RoughImage1;
        RoughCalibrated1(k) = RoughCalibrated1(k).*squeeze(CalibrationFactor(1,i,:));
        subplot(4,4,2*(i-1)+1);imagesc(RoughCalibrated1);title(title1);
        RoughCalibrated2 = RoughImage2;
        RoughCalibrated2(k) = RoughCalibrated2(k).*squeeze(CalibrationFactor(1,i,:));
        subplot(4,4,2*(i-1)+2);imagesc(RoughCalibrated2);title(title2);
    end
    
%     figure;imagesc(RoughImage1);
%     RoughImage1(k) = RoughImage1(k).*squeeze(CalibrationFactor(1,i,:));
%     figure;imagesc(RoughImage1);

    VectorImages(1,i,:) = RoughImage1(k);
    VectorImages(2,i,:) = RoughImage2(k);
%     VectorImages(1,i,:) = RoughImage1(k).*squeeze(CalibrationFactor(1,i,:));
%     VectorImages(2,i,:) = RoughImage2(k).*squeeze(CalibrationFactor(2,i,:));
end

end

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