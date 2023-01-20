function [VectorImage1,VectorImage2] = get_SXRImage(date,number,SXRfilename,filter)
% % 画像を切り取る
N_projection = 80;
% VectorImages = CutImage(date,shot,N_projection/230,false);
CheckFlag = false;
VectorImages = CutImage(date,N_projection/240,CheckFlag,SXRfilename,filter);
VectorImages1 = squeeze(VectorImages(1,:,:));
VectorImages2 = squeeze(VectorImages(2,:,:));
VectorImage1 = VectorImages1(number,:);
VectorImage2 = VectorImages2(number,:);
end

function VectorImages = CutImage(date,resolution,check_flag,SXRfilename,filter)


% pathname = '/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/SXR_Images/';
% foldername = strcat(pathname,num2str(date));
% shotname = strcat('/shot',num2str(shotnumber),'.tif');
% 
% if date <= 210924
%     shotname = strcat('/shot',num2str(shotnumber),'.tif');
% else
%     if  shotnumber < 10
%         shotname = strcat('/',num2str(date),'00',num2str(shotnumber),'.tif');
%     elseif shotnumber < 100
%         shotname = strcat('/',num2str(date),'0',num2str(shotnumber),'.tif');
%     elseif shotnumber < 1000
%         shotname = strcat('/',num2str(date),num2str(shotnumber),'.tif');
%     else
%         disp('More than 999 shots! You need some rest!!!')
%         return
%     end
% end

% filename = strcat(foldername,shotname);

RawImage = imread(SXRfilename);

if filter
    figure;imagesc(RawImage);
    [RawImage,~] = imnlmfilt(RawImage,'SearchWindowSize',91,'ComparisonWindowSize',15);
    figure;imagesc(RawImage);
% else
%     RawImage = RawImage;
end
% whos RawImage

% IW = 115; %Image Width
IW = 120; %Image Width

TimeRapsImage = zeros(2,8,2*IW,2*IW);

Center = zeros(2,8,2);
% 
% % set the centers of each image
% Center(2,1,:) = [500,1585];
% Center(2,2,:) = [1035,1620];
% Center(2,3,:) = [1055,1245];
% Center(2,4,:) = [505,1225];
% Center(2,5,:) = [515,840];
% Center(2,6,:) = [1060,860];
% Center(2,7,:) = [1060,485];
% Center(2,8,:) = [515,455];
% 
% Center(2,:,1) = Center(2,:,1) + 25;
% if date >= 210310
%     Center(2,:,2) = Center(2,:,2) - 10;
% elseif date == 210309
%     Center(2,:,2) = Center(2,:,2) - 20;
% end
% 
% Center(1,:,1) = Center(2,:,1) - 240;
% Center(1,:,2) = Center(2,:,2) - 15;

if date < 210902
    Center(2,1,:) = [500,1585];
    Center(2,2,:) = [1035,1620];
    Center(2,3,:) = [1055,1245];
    Center(2,4,:) = [505,1225];
    Center(2,5,:) = [515,840];
    Center(2,6,:) = [1060,860];
    Center(2,7,:) = [1060,485];
    Center(2,8,:) = [515,455];
    Center(2,:,1) = Center(2,:,1) + 25;
    if date >= 210310
        Center(2,:,2) = Center(2,:,2) - 10;
    elseif date == 210309
        Center(2,:,2) = Center(2,:,2) - 20;
    end
    Center(1,:,1) = Center(2,:,1) - 240;
    Center(1,:,2) = Center(2,:,2) - 15;
elseif (date >= 210902) && (date < 210916)
    Center(2,1,:) = [480,1585];
    Center(2,2,:) = [1015,1620];
    Center(2,3,:) = [1030,1245];
    Center(2,4,:) = [480,1225];
    Center(2,5,:) = [490,840];
    Center(2,6,:) = [1040,870];
    Center(2,7,:) = [1040,480];
    Center(2,8,:) = [500,450];
    Center(2,:,1) = Center(2,:,1) + 25;
    Center(1,:,1) = Center(2,:,1) - 235;
    Center(1,:,2) = Center(2,:,2) - 10;
elseif (date >= 210916) && (date < 210924)
    Center(2,1,:) = [490,1585];
    Center(2,2,:) = [1030,1610];
    Center(2,3,:) = [1040,1240];
    Center(2,4,:) = [490,1220];
    Center(2,5,:) = [500,840];
    Center(2,6,:) = [1045,870];
    Center(2,7,:) = [1055,480];
    Center(2,8,:) = [510,450];
    Center(1,:,1) = Center(2,:,1) - 235;
    Center(1,:,2) = Center(2,:,2) - 10;
elseif (date >= 210924) && (date < 211223)
    Center(2,1,:) = [510,1585];
    Center(2,2,:) = [1060,1610];
    Center(2,3,:) = [1060,1240];
    Center(2,4,:) = [520,1220];
    Center(2,5,:) = [520,840];
    Center(2,6,:) = [1060,870];
    Center(2,7,:) = [1070,480];
    Center(2,8,:) = [530,450];
    Center(1,:,1) = Center(2,:,1) - 235;
    Center(1,:,2) = Center(2,:,2) - 10;
%     Center(2,1,:) = [490,1585];
%     Center(2,2,:) = [1030,1610];
%     Center(2,3,:) = [1040,1240];
%     Center(2,4,:) = [490,1220];
%     Center(2,5,:) = [500,840];
%     Center(2,6,:) = [1045,870];
%     Center(2,7,:) = [1055,480];
%     Center(2,8,:) = [510,450];
%     Center(1,:,1) = Center(2,:,1) - 235;
%     Center(1,:,2) = Center(2,:,2) - 10;
elseif date >= 211223
    Center(2,1,:) = [510,1565];
    Center(2,2,:) = [1040,1590];
    Center(2,3,:) = [1060,1220];
    Center(2,4,:) = [510,1200];
    Center(2,5,:) = [520,815];
    Center(2,6,:) = [1070,855];
    Center(2,7,:) = [1075,460];
    Center(2,8,:) = [535,430];
    Center(1,:,1) = Center(2,:,1) - 235;
    Center(1,:,2) = Center(2,:,2) - 10;
end

% disp(Center);

for i = 1:8
    UpRangeV = Center(1,i,1)-IW+1:Center(1,i,1)+IW;
    UpRangeH = Center(1,i,2)-IW+1:Center(1,i,2)+IW;
    DownRangeV = Center(2,i,1)-IW+1:Center(2,i,1)+IW;
    DownRangeH = Center(2,i,2)-IW+1:Center(2,i,2)+IW;
    TimeRapsImage(1,i,:,:) = RawImage(UpRangeV,UpRangeH,1);
    TimeRapsImage(2,i,:,:) = RawImage(DownRangeV,DownRangeH,1);
    if date >= 211223
        TimeRapsImage(1,i,:,:) = rot90(TimeRapsImage(1,i,:,:),2);
        TimeRapsImage(2,i,:,:) = rot90(TimeRapsImage(2,i,:,:),2);
    end
end

BackGround = cast(RawImage(1050-IW+1:1050+IW,120-IW+1:120+IW,1),'double');

k = FindCircle(IW*resolution);
VectorImages = zeros(2,8,numel(k));

if check_flag
    f=figure;
%     f.WindowState = 'maximized';
    f.Position = [900,200,500,500];
end

for i=1:8
    SingleImage1 = squeeze(TimeRapsImage(1,i,:,:))-BackGround;
    SingleImage2 = squeeze(TimeRapsImage(2,i,:,:))-BackGround;
%     SingleImage1 = squeeze(TimeRapsImage(1,i,:,:));
%     SingleImage2 = squeeze(TimeRapsImage(2,i,:,:));
%     GrayImage1 = flipud(SingleImage1(:,:)); %反転なし
%     GrayImage2 = flipud(SingleImage2(:,:));
    GrayImage1 = fliplr(SingleImage1(:,:)); %反転なし
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
        figure(f);
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