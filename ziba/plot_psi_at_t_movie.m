function [] = plot_psi_at_t_movie(B_z,r_probe,z_probe,time,fitting,fill,fixed_Clayer,show_probe)
% plot psi in rz plane
% input:
%   3d array of double: B_z (r,z,t), offsetted at zero and smoothed
%   1d array of double: r_probe, locations of probes along r
%   1d array of double: z_probe, locations of probes along z
%   integer: t, time of interest (us)
%   boolean: fitting, option for enabling cubic spline fitting 
%   to fit the the 8x28 mesh to 50x50 mesh
%   boolean: fixed_Clayer, option for enabling 
%            true -> fixed value between each contour layer
%            false-> fixed number of total layers = layer_num
%   boolean: show_probe, option for showing location of pickup coils

%data画像保存先/video保存名
pathname.save='C:\Users\uswk0\OneDrive\デスクトップ\data\movie_out\'; %保存先
video_name = 'out';

layer_num = 100;
% a uniform mesh for interpolation
[probe_mesh_z,probe_mesh_r] = meshgrid(z_probe,r_probe);

if fitting
    z_space = linspace(z_probe(1),z_probe(end),50);
    r_space = linspace(r_probe(1),r_probe(end),50);
    [psi_mesh_z,psi_mesh_r] = meshgrid(z_space,r_space);
else
    [psi_mesh_z,psi_mesh_r] = meshgrid(z_probe,r_probe);
end

if fixed_Clayer
    max_psi_color = 0.0400;
    min_psi_color = -0.0200;
    layer_resolution = 0.0002;
    contour_layer =  min_psi_color:layer_resolution:max_psi_color;
else
    contour_layer = layer_num;
end

for i = time
figure('Position', [0 0 1200 600])
hold on

psi = get_psi(B_z,r_probe,i);

if fitting
    psi = griddata(z_probe,r_probe,psi,psi_mesh_z,psi_mesh_r,'cubic');
end
if fill
    contourf(psi_mesh_z,psi_mesh_r,psi,contour_layer,'Fill','on');
else
    contourf(psi_mesh_z,psi_mesh_r,psi,contour_layer,'--b','Fill','off');
end

if show_probe
    plot(probe_mesh_z,probe_mesh_r,'*','color','k','markersize',3)
end
    
title(strcat(num2str(i),' us','; Psi (Wb)'));
xlabel('z (m)');
ylabel('r (m)');
ax = gca;
ax.FontSize = 18; 


hold off

%画像の作成
saveas(gcf,strcat(pathname.save,'\time',num2str(i),'_shot',num2str(date),'.jpg'))
close


daspect([1 1 1])

end



%%%%%%%%%%%%　動画作成用コード %%%%%%%%%%%%%%
workingDir = pathname.save;
imageNames = dir(fullfile(workingDir,'*.jpg'));
imageNames = {imageNames.name}';

outputVideo = VideoWriter(fullfile(workingDir,video_name),"MPEG-4");
%再生速度の変更(defaultの30は早すぎ)
outputVideo.FrameRate = 5;

open(outputVideo)
%writeAnimation(outputVideo)

%イメージ シーケンス内で繰り返し、各イメージを読み取り、それをビデオに書き込みます。

for ii = 1:length(imageNames)
   img = imread(fullfile(workingDir,imageNames{ii}));
   writeVideo(outputVideo,img)
end
%ビデオ ファイルを完成します。

close(outputVideo)
%最終的なビデオの表示
%リーダー オブジェクトを作成します。
shuttleMp4 = VideoReader(fullfile(workingDir,strcat(video_name,'.mp4')));

%ビデオ フレームから、MATLAB® ムービー struct を作成します。
ii = 1;
while hasFrame(shuttleMp4)
   mov(ii) = im2frame(readFrame(shuttleMp4));
   ii = ii+1;
end

%ビデオの幅および高さに基づいて現在の Figure および Axes のサイズを調整し、ムービーの最初のフレームを表示します。
figure 
imshow(mov(1).cdata, 'Border', 'tight')

%ビデオのフレーム レートでムービーを 1 回再生します。
movie(mov,1,shuttleMp4.FrameRate)

