clear
%%%%%%%%%%%%%%%%%%%%%%%%
%  288chDopplerとmrdの磁気面を重ねてプロットして保存するコード
%　288chDopplerはIDLで作ったsavのデータをあらかじめI:\makimitsu\yyddmmに保存してあるものを読み込む
%
%%%%%%%%%%%%%%%%%%%%%%%%

yourname = 'C:\Users\Moe Akimitsu\';
f = fullfile(yourname,'Documents','GitHub','test-open');
addpath(genpath(f));
%addpath(fullfile(yourname,'Documents','GitHub','test-open','pcb'));
%f = fullfile(yourname,'Documents','GitHub','SXR_test');
%addpath(f);
%%%適宜変更


DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';
T=getTS6log(DOCID);% ログのテーブルを取得
 
shotlist=[2692:2950];
subT=T(shotlist,:);
IDXlist=shotlist(isfinite(subT.DopplerDelay)&isfinite(subT.d_tacq));
%IDX=IDXlist(1,88);
for IDX=IDXlist%(1,49:end)  %
date=T.date(IDX);
shot=T.shot(IDX);
TF_shot=T.TFoffset(IDX);
if isnan(T.EF_A_(IDX))%%NaNでないことを確認（ログが空白だとNaNになる）
    i_EF=150;
else  %NaNなら150をとりあえず代入、記入されているときはその値を使う
    i_EF=T.EF_A_(IDX);
end
t=T.DopplerDelay(IDX);
n=50;
[B_z,r_probe,z_probe,ch_dist,data,data_raw,shot_num] = get_B_z(date,TF_shot,shot,offset_TF,i_EF,folder_path);

B_z = B_z([2,3,4,6,7,8],2:end,:);
data = data([2,3,4,6,7,8],2:end,:);
z_probe = z_probe(2:end);
ch_dist = ch_dist([2,3,4,6,7,8],2:end);
r_probe = r_probe([2,3,4,6,7,8]);

%%ファイルへのパスを作る
%それぞれのPCから共有フォルダまでのパスはそれぞれ異なるので各自で設定
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス
pathname.fourier='I:';%md0までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath

%共有フォルダ以下から目的ショットのファイルを探す
%filepath.rgw=strcat(pathname.ts3u, '\', string(date),'\' ...
%    ,string(date),num2str(shot,'%03i'),'.rgw');
%filepath.D288=dir(strcat(pathname.fourier,'\Doppler\288CH\20',string(date),'\*shot',num2str(shot),'*.asc'));
if shot<10
    filepath.D288=dir(strcat(pathname.fourier,'\makimitsu\',string(date),'\doppler2D_shot',num2str(shot),'_.sav'));
else
    filepath.D288=dir(strcat(pathname.fourier,'\makimitsu\',string(date),'\doppler2D_shot',num2str(shot),'.sav'));
end    
    %filepath.Dhighspeed=dir(strcat(pathname.NIFS,'\Doppler\Photron\',string(date),'\**\*shot',num2str(shot),'*.tif'));
    %filepath.SXR=strcat(pathname.NIFS,'\X-ray\',string(date),'\shots\',string(date),num2str(shot,'%03i'),'.tif');
if numel(filepath.D288)==0
    continue
end

if isfile( fullfile(filepath.D288.folder,filepath.D288.name))
    psi = get_psi(B_z,r_probe,t);
    z_space = linspace(z_probe(1),z_probe(end),50);
    r_space = linspace(r_probe(1),r_probe(end),50);
    [psi_mesh_z,psi_mesh_r] = meshgrid(z_space,r_space);
    [probe_mesh_z,probe_mesh_r] = meshgrid(z_probe,r_probe);
    psi = griddata(z_probe,r_probe,psi,psi_mesh_z,psi_mesh_r,'cubic');
    restore_idl( fullfile(filepath.D288.folder,filepath.D288.name),'lowercase','create'); %.savファイルを読み込む
    
    f = figure;
    % f.WindowState = 'maximized';
    f.Units = 'normalized';
    f.Position = [0.1,0.2,0.8,0.6];
    pos1 = [0.07,0.2,0.35,0.6];
    pos2 = [0.58,0.2,0.35,0.6];

    subplot('Position',pos1);
    contourf(doppler.z,doppler.yy,doppler.emission,10,'LineStyle','none')
    colormap(jet)
    axis image
    axis tight manual
    colorbar('Location','eastoutside')
    hold on
    contour(psi_mesh_z,psi_mesh_r,psi,30,'-b');
    hold off
    title(string(t)+'us,emiision')
    xlabel('z')
    ylabel('r')
    caxis([-2e5,2e5])

    subplot('Position',pos2);
    contourf(doppler.z,doppler.yy,doppler.ti_2d,30,'LineStyle','none')
    colormap(jet)
    axis image
    axis tight manual
    colorbar('Location','eastoutside')
    hold on
    contour(psi_mesh_z,psi_mesh_r,psi,30,'-b');
    hold off
    caxis([-300,300])
    title(string(t)+'us,ti')
    xlabel('z')
    ylabel('r')
    filename = strcat('I:\makimitsu\',num2str(date),'\Doppler_',num2str(date),num2str(shot,'%03i'),'_',num2str(t),'us');
    saveas(gcf,strcat(filename,'.png'))
    close
end


% %[B_z,r_probe,z_probe,ch_dist,B_z_return,data_return,shot_num] 
% [B_z,r_probe,z_probe,ch_dist,data,data_raw,shot_num]= get_B_z(date,TF_shot,shot,offset_TF,T.EF_A_(IDX),pathname.ts3u);
% B_z = B_z([2,3,4,6,7,8],2:end,:);
% data = data([2,3,4,6,7,8],2:end,:);
% z_probe = z_probe(2:end);
% ch_dist = ch_dist([2,3,4,6,7,8],2:end);
% r_probe = r_probe([2,3,4,6,7,8]);



%plot_B_z_in_time(B_z,ch_dist,350,600);
%plot_psi_SXR_multi(B_z,r_probe,z_probe,date,shot,layer,area,start,exposure,SXRfilename)
%plot_psi_SXR_multi(B_z,r_probe,z_probe,date,shot,true,true,T.Period_StartTime_(IDX),2,filepath.SXR)
end
%

