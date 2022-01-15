clear
%%%%%%%%%%%%%%%%%%%%%%%%
%  SXRとpcbデータを重ねてプロットして保存するコード
%　 test-openで実行する。ログから読み込んで、それぞれのパスを設定するコード
%%　参考　https://jp.mathworks.com/help/matlab/matlab_prog/access-data-in-a-table.html
%%%%%%%%%%%%%%%%%%%%%%%%

yourname = 'C:\Users\Moe Akimitsu\';
f = fullfile(yourname,'Documents','GitHub','test-open');
addpath(genpath(f));
%addpath(fullfile(yourname,'Documents','GitHub','test-open','pcb'));
f = fullfile(yourname,'Documents','GitHub','SXR_test');
addpath(f);
%%%適宜変更


DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';
T=getTS6log(DOCID);% ログのテーブルを取得

shotlist=[2692:2950];
subT=T(shotlist,:);
IDXlist=shotlist(isfinite(subT.Period_StartTime_)&isfinite(subT.d_tacq));
%IDX=IDXlist(1,88);
for IDX=IDXlist%(1,1:5)
date=T.date(IDX);
shot=T.shot(IDX);
TF_shot=T.TFoffset(IDX);
offset_TF=true;
start=T.Period_StartTime_(IDX);
if isnan(T.EF_A_(IDX))
    i_EF=150;
else  
    i_EF=T.EF_A_(IDX);
end



%%ファイルへのパスを作る
%それぞれのPCから共有フォルダまでのパスはそれぞれ異なるので各自で設定
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス
pathname.fourier='I:';%md0までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath

%共有フォルダ以下から目的ショットのファイルを探す
filepath.rgw=strcat(pathname.ts3u, '\', string(date),'\' ...
    ,string(date),num2str(shot,'%03i'),'.rgw');
%filepath.D288=dir(strcat(pathname.fourier,'\Doppler\288CH\20',string(date),'\*shot',num2str(shot),'*.asc'));
%filepath.Dhighspeed=dir(strcat(pathname.NIFS,'\Doppler\Photron\',string(date),'\**\*shot',num2str(shot),'*.tif'));
%filepath.SXR=strcat(pathname.NIFS,'\X-ray\',string(date),'\shots\',string(date),num2str(shot,'%03i'),'.tif');
times = start:5:(start+5*7);
cd 'I:\makimitsu\211223';
for t = times
        
        filename = strcat('SXR_',num2str(date),num2str(shot,'%03i'),'_',num2str(t),'us');
        filemat=strcat(filename,'.mat');

if isfile(filemat)

   
     load(filemat)
    %[B_z,r_probe,z_probe,ch_dist,B_z_return,data_return,shot_num] 
    [B_z,r_probe,z_probe,ch_dist,data,data_raw,shot_num]= get_B_z(date,TF_shot,shot,offset_TF,T.EF_A_(IDX),pathname.ts3u);
    B_z = B_z([2,3,4,6,7,8],2:end,:);
    data = data([2,3,4,6,7,8],2:end,:);
    z_probe = z_probe(2:end);
    ch_dist = ch_dist([2,3,4,6,7,8],2:end);
    r_probe = r_probe([2,3,4,6,7,8]);
    z_space = linspace(z_probe(1),z_probe(end),50);
    r_space = linspace(r_probe(1),r_probe(end),50);
    [psi_mesh_z,psi_mesh_r] = meshgrid(z_space,r_space);
    
    

      
    psi = get_psi(B_z,r_probe,t);
    psi = griddata(z_probe,r_probe,psi,psi_mesh_z,psi_mesh_r,'cubic');
    
    f = figure;
    % f.WindowState = 'maximized';
    f.Units = 'normalized';
    f.Position = [0.1,0.2,0.8,0.6];
    pos1 = [0.07,0.2,0.35,0.6];
    pos2 = [0.58,0.2,0.35,0.6];
    
    % subplot(2,1,1);
    subplot('Position',pos1);
    [~,h1] = contourf(SXR_mesh_z1,SXR_mesh_r1,EE1,40);
    h1.LineStyle = 'none';
    % caxis([0,0.5]);
    % caxis([0,0.2]);
    c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
    hold on
    [~,hp1]=contourf(psi_mesh_z,psi_mesh_r,psi,30,'-k','Fill','off');

    title(strcat(num2str(t),' us'));
    xlabel('z [m]');
    ylabel('r [m]');
    ax = gca;
    ax.FontSize = 18; 
    hold off
        axis image
        axis tight manual
    % subplot(2,1,2);
    subplot('Position',pos2);

    [~,h2] = contourf(SXR_mesh_z2,SXR_mesh_r2,EE2,40);
    h2.LineStyle = 'none';
    % caxis([0,0.15]);
    % caxis([0,0.2]);
    c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
    hold on
    [~,hp2]=contourf(psi_mesh_z,psi_mesh_r,psi,30,'-k','Fill','off');
    title(strcat(num2str(t),' us'));
    xlabel('z [m]');
    ylabel('r [m]');
    ax = gca;
    ax.FontSize = 18; 
        axis image
        axis tight manual
    hold off
    %plot_B_z_in_time(B_z,ch_dist,350,600);
    %plot_psi_SXR_multi(B_z,r_probe,z_probe,date,shot,layer,area,start,exposure,SXRfilename)
    %plot_psi_SXR_multi(B_z,r_probe,z_probe,date,shot,true,true,T.Period_StartTime_(IDX),5,true,filepath.SXR)
    saveas(gcf,strcat(filename,'_all.png'))
    close
end
end
end
%DopplerDelay=T.DopplerDelay;

