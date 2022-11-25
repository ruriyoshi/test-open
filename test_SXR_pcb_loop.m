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

%shotlist=[2692:2950];
shotlist=[2911:2950];
subT=T(shotlist,:);
IDXlist=shotlist(isfinite(subT.Period_StartTime_)&isfinite(subT.d_tacq));
%IDX=IDXlist(1,88);
for IDX=IDXlist(1,6:end)
date=T.date(IDX);
shot=T.shot(IDX);
d_tacq=T.d_tacq(IDX);
d_tacqTF=T.TFdtacq(IDX);
if isnan(T.EF_A_(IDX))
    i_EF=150;
else  
    i_EF=T.EF_A_(IDX);
end
%disp(['EF=', num2str(T.EF_A_(IDX))])
%%NaNでないことを確認（ログが空白だとNaNになる）
trange=T.Period_StartTime_(IDX):5:T.Period_StartTime_(IDX)+7*5;
n=50;
[grid2D, data2D] = pcbdata(date, d_tacq,d_tacqTF,trange, [], n,i_EF);
if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ
    continue
end
%     maxrange=max(abs(data2D.Jt),[],'all');
% [psimid,mid]=min(data2D.psi,[],2);
% [opoint,p]=islocalmin(psimid,1);
% [xpoint,~]=islocalmax(psimid,1);
% [xp_psi,maxxp]=max(squeeze(psimid),[],1);
% onum=squeeze(sum(opoint,1));
% trange(onum~=0)
    %maxrange=2e6;

%%磁気面プロット
% figure
% start=8; %460+?
% for m=1:10 
%     i=start+m;
%     t=trange(i);
%     subplot(2,5,m)
%     contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Jt(:,:,i),10,'LineStyle','none')
%     colormap(jet)
%     axis image
%     axis tight manual
%     %     xlim([-0.02 0.02])
%     %     ylim([0.12 0.27])
%     caxis([-maxrange,maxrange])
%     colorbar('Location','eastoutside')
%     %zlim([-1 1])
%     %colormap(bone)
%     hold on
%     plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
%     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),50,'black')
%     plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"ro")
%     plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"rx")
%     hold off
%     title(string(t)+'us')
%     xlabel('z')
%     ylabel('r')
% end


%%ファイルへのパスを作る
%それぞれのPCから共有フォルダまでのパスはそれぞれ異なるので各自で設定
%pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス
%pathname.fourier='I:';%md0までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath

%共有フォルダ以下から目的ショットのファイルを探す
%filepath.rgw=strcat(pathname.ts3u, '\', string(date),'\' ...
%    ,string(date),num2str(shot,'%03i'),'.rgw');
%filepath.D288=dir(strcat(pathname.fourier,'\Doppler\288CH\20',string(date),'\*shot',num2str(shot),'*.asc'));
%filepath.Dhighspeed=dir(strcat(pathname.NIFS,'\Doppler\Photron\',string(date),'\**\*shot',num2str(shot),'*.tif'));
filepath.SXR=strcat(pathname.NIFS,'\X-ray\',string(date),'\shots\',string(date),num2str(shot,'%03i'),'.tif');

if isfile(filepath.SXR)
    %plot_pcb_SXR_multi(data2D.psi,grid2D.rq,grid2D.zq,date,shot,true,false,T.Period_StartTime_(IDX),5,true,filepath.SXR)
plot_pcb_SXR_multi(data2D.psi,grid2D.rq,grid2D.zq,date,shot,true,false,T.Period_StartTime_(IDX),5,true,filepath.SXR)
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
%DopplerDelay=T.DopplerDelay;

