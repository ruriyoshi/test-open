function [VectorImage1,VectorImage2] = get_SXRImage(Date,number,SXRImageFilePath,FileterProcessFlag)
% この関数は，撮影日時(date)からファイバー中心を確定し，X線画像(SXRImageFilePath)を切り取り，再構成直前のベクトル画像を返します．

% date = 230119;
% number = 5;
% SXRfilename = '/Users/shinjirotakeda/Library/CloudStorage/OneDrive-TheUniversityofTokyo/Documents/SXR_Images/230119/shot024.tif';
% filter = false;

% N_Projectionは，再構成直前のベクトル画像を二次元に変換したときの一辺の長さに相当します．つまり，計算の負荷に直結するのがN_Projectionです．
% 一方RawImageFiberHalfWidthは，RawImageにおけるファイバー円の半径に相当します．これはファイバーバンドルと高速度カメラの光学的な位置関係によって決定されます．
N_Projection = 80;
% VectorImages = CutImage(date,shot,N_projection/230,false);
CheckFlag = true;
% CheckFlag = false;
VectorImages = CutImage(Date,N_Projection,CheckFlag,SXRImageFilePath,FileterProcessFlag);
VectorImages1 = squeeze(VectorImages(1,:,:));
VectorImages2 = squeeze(VectorImages(2,:,:));
VectorImage1 = VectorImages1(number,:);
VectorImage2 = VectorImages2(number,:);
end

function VectorImages = CutImage(Date,N_Projection,check_flag,SXRImageFilePath,FilterProcessFlag)

RawImage = imread(SXRImageFilePath);

if FilterProcessFlag
    figure;imagesc(RawImage);
    [RawImage,~] = imnlmfilt(RawImage,'SearchWindowSize',91,'ComparisonWindowSize',15);
    figure;imagesc(RawImage);
end

FiberCalibrationImagePath = strcat('G:/My Drive/X-ray/Data/TIF/',num2str(Date),'/PositionCheck.tif');
CalibrationImage = imread(FiberCalibrationImagePath);
% 以降，if(date)分岐によって既にファイバー中心が判明しているものはそれを代入，判明していないものはFindFiber関数で中心を探します．
% centersは(16,2)行列で，例えばcenters(5,:)はRawImageがもつ16個の円画像のうち右上から蛇順に数えて5番目の円の中心の(row,column)を表します．
% 理想的にはIW=radiiですが，ここではradiiを信用せずファイバーの直径IWを手動で設定しています．
if Date == 230315
    centers = 1.0e3*[...
        0.2815,1.5525;0.4834,1.5468;0.8174,1.5699;1.0084,1.5579;...
        1.0215,1.1871;0.8132,1.1972;0.4748,1.1887;0.2770,1.1942;...
        0.2766,0.8066;0.4750,0.7979;0.8203,0.8339;1.0226,0.8280;...
        1.0169,0.4307;0.8187,0.4377;0.4819,0.4108;0.2794,0.4245];
    RawImageFiberHalfWidth = 95;
elseif Date == 230119
    centers = 1.0e3*[...
        0.2815,1.5525;0.5234,1.5168;0.8174,1.5699;1.0584,1.5279;...
        1.0665,1.1621;0.8532,1.1972;0.5248,1.1487;0.3370,1.1542;...
        0.3366,0.7766;0.5250,0.7679;0.8753,0.7939;1.0676,0.7830;...
        1.0669,0.4107;0.8787,0.4177;0.5319,0.3808;0.3394,0.3945];
    RawImageFiberHalfWidth = 90;
elseif Date == 230316
    [centers,radii]=FindFibers(CalibrationImage,[90,130]);
    RawImageFiberHalfWidth = 90;
    centers(:,1) = centers(:,1) + 40;
    centers(:,2) = centers(:,2) - 20;
else
    [centers,radii]=FindFibers(CalibrationImage,[70,130]);
    RawImageFiberHalfWidth = 85;
end
% 上にずらすときは1つ目の数字を増やす
% 右にずらすときは2つめの数字を増やす

% IWをradiiから求めようとする場合は，次の行を使用します．
% IW = round(mean(radii));
centers = round(centers);
Center = zeros(2,8,2);
% TimeRapsImageはファイバ画像を二次元で格納する変数で，例えばTimeRapsImage(1,5,:,:)は一枚目のX線フィルタ，時系列で5番目のファイバ二次元画像を表します．
TimeRapsImage = zeros(2,8,2*RawImageFiberHalfWidth,2*RawImageFiberHalfWidth);
% Centerは(2 8 2)行列で，例えばCenter(1,5,:)は一枚目のX線フィルタ，時系列で5番目のファイバ中心の座標(row,column)を示します．
% centersには16枚の円の中心が蛇順に格納されているので，これをX線フィルタごと，時系列順に整理します．
Center(1,:,:) = centers([1,3,6,8,9,11,14,16],:);
Center(2,:,:) = centers([2,4,5,7,10,12,13,15],:);
% Center(2,1,:) = [530,1520];
% Center(2,2,:) = [1060,1530];
% Center(2,3,:) = [1070,1160];
% Center(2,4,:) = [530,1150];
% Center(2,5,:) = [530,770];
% Center(2,6,:) = [1070,780];
% Center(2,7,:) = [1070,410];
% Center(2,8,:) = [530,380];
% Center(1,:,1) = Center(2,:,1) - 190;
% Center(1,:,2) = Center(2,:,2) + 5;

% disp(Center);
% whos RawImage
% whos TimeRapsImage

% RawImageを矩形に切り取り，TimeRapsImageに格納します．
for i = 1:8
    UpFiberImageWidth_row = Center(1,i,1)-RawImageFiberHalfWidth+1:Center(1,i,1)+RawImageFiberHalfWidth;
    UpFiberImageWidth_column = Center(1,i,2)-RawImageFiberHalfWidth+1:Center(1,i,2)+RawImageFiberHalfWidth;
    DownFiberImageWidth_row = Center(2,i,1)-RawImageFiberHalfWidth+1:Center(2,i,1)+RawImageFiberHalfWidth;
    DownFiberImageWidth_column = Center(2,i,2)-RawImageFiberHalfWidth+1:Center(2,i,2)+RawImageFiberHalfWidth;
    TimeRapsImage(1,i,:,:) = RawImage(UpFiberImageWidth_row,UpFiberImageWidth_column,1);
    TimeRapsImage(2,i,:,:) = RawImage(DownFiberImageWidth_row,DownFiberImageWidth_column,1);
    TimeRapsImage(1,i,:,:) = rot90(TimeRapsImage(1,i,:,:),2);
    TimeRapsImage(2,i,:,:) = rot90(TimeRapsImage(2,i,:,:),2);
end
% RawImageのうちファイバが投影されていない部分(つまり理想的には0で合ってほしい部分)を背景ノイズとして扱います．
BackGround = cast(RawImage(1050-RawImageFiberHalfWidth+1:1050+RawImageFiberHalfWidth,120-RawImageFiberHalfWidth+1:120+RawImageFiberHalfWidth,1),'double');
BackGround = ones(size(BackGround))*mean(BackGround,'all');
% TimeRapsImageの中からファイバが投影された円だけを取り出し，画質を変換し，二次元を一次元に変換し，VectorImagesに格納します．
% kはN_Projectionを一辺にもつ正方形の中で，内接円の内側に位置するピクセルのインデックスを持ちます．
k = FindCircle(N_Projection/2);
VectorImages = zeros(2,8,numel(k));
if check_flag
    f1=figure;
%     f.WindowState = 'maximized';
    f1.Position = [900,200,500,500];
%     f2=figure;
%     f.WindowState = 'maximized';
%     f2.Position = [900,200,500,500];
end

% 関数clc_CalibrationFactorはファイバ時系列全ての強度の平均値に対する各ピクセルの強度を計算し保存します．これは放電間のX線強度を正規化する為に行いますが，現在は使用しません．
CalibrationPath = strcat('G:/My Drive/X-ray/Data/TIF/',num2str(Date),'/CalibrationFactor.mat');
if exist(CalibrationPath,'file')
    load(CalibrationPath,'CalibrationFactor');
else
    CalibrationFactor = clc_CalibrationFactor(Date,N_Projection);
end

% RawImageから切り取られた円の直径ピクセルはRawImageFiberHalfWidthです．これを直径ピクセルN_Projectionに変換します．
% SingleImage,GrayImage,RoughImageは全て二次元画像です．
resolution = N_Projection/(RawImageFiberHalfWidth*2);
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
%         figure(f2);
%         RoughCalibrated1 = RoughImage1;
%         RoughCalibrated1(k) = RoughCalibrated1(k).*squeeze(CalibrationFactor(1,i,:));
%         subplot(4,4,2*(i-1)+1);imagesc(RoughCalibrated1);title(title1);
%         RoughCalibrated2 = RoughImage2;
%         RoughCalibrated2(k) = RoughCalibrated2(k).*squeeze(CalibrationFactor(1,i,:));
%         subplot(4,4,2*(i-1)+2);imagesc(RoughCalibrated2);title(title2);
    end
%     figure;imagesc(RoughImage1);
%     RoughImage1(k) = RoughImage1(k).*squeeze(CalibrationFactor(1,i,:));
%     figure;imagesc(RoughImage1);
    % ここでRoughImageという二次元画像がVectorImageという一次元画像に変換されます．
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