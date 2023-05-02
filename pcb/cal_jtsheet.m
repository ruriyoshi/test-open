%%%%%%%%%%%%%%%%%%%%%%%%
%TF vs Vin, Voutの解析
%%%%%%%%%%%%%%%%%%%%%%%%

%%ファイルへのパスを作る
%それぞれのPCから共有フォルダまでのパスはそれぞれ異なるので各自で設定

%%%%(1)spread sheetから ログのテーブルを取得してTに格納
%Github/test-open/getTS6log.mを使用
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);

%%%%%ここが各PCのパス
%環境変数を設定していない場合はパスを''内に全て記入する（使用しないパスは空白''で良い）
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.save=getenv('save_path'); %保存先

%%%%(2)ログから解析したいデータを検索
%Github/test-open/searchlog.mを使用

% node='date';  % 【input】検索する列の名前. T.Properties.VariableNamesで一覧表示できる
%  pat=211223;   % 【input】検索パターン（数値なら一致検索、文字なら含む検索）　
% 
% searchlog(T,node,pat); % ログのテーブルから当てはまるものを抽出した新しいテーブルを作成

%%%%(3)指定したshotの解析
IDXlist=[2911:2913 2925 2926 2927 2931 2933 2947:2950 2942 2943 2946];

% shot=struct('a',2911:2913,'b',[2925 2926 2927],'c',[2931 2933],'d',2947:2950,'e',[2942 2943 2946]);
% IDXlist=shot.a; %【input】テーブルから解析したいshot番号を抽出して入力
%(a)211223/IDX2870:2920(42:44)/2911:2913/dtacq3592:3594
%(b)211224/IDX2922:2950(4,5,8)/2925,2926,2927/dtacq3606,3607,3609
%(c)211224/IDX2922:2950(10,12)/2931,2933/dtacq3611,3613
%(d)211224/IDX2922:2950(26:29)/2947:2950/dtacq3626-3629
%(e)211224/IDX2922:2950(21,22,25)/2942,2943,2946/dtacq3621,3622,3625

% tcal=[472 478; 474 479; 474 481; 
%     476 484; 476 482; 476 483; 
%     485 499; 485 499;
%     475 483; 474 483; 473 483; 479 489;
%     485 496; 484 494; 483 495]; %15*2

% figure
for i=4 %1:numel(IDXlist)
    IDX=IDXlist(i);
    [xEt,xJt,eta,trange]=def_sheet(T, IDX,pathname);
    %plot(trange(1:end-1),Et_midmax)
    %plot(trange,xEt,'--o')
    %plot(trange,xJt,'-')
    %hold on
end
%hold off
%title(strcat('IDX=',num2str(IDXlist(42)),':',num2str(IDXlist(44))))
%xlim([trange(1) trange(end-1)])
%xlabel('Time [us]')
%ylabel('Jt at Xpoint [A/m^2]')
%ylabel('Et at Xpoint [V/m]')
%legend(strcat('IDX=',num2str(IDXlist(1))),strcat('IDX=',num2str(IDXlist(2))),strcat('IDX=',num2str(IDXlist(3))),strcat('IDX=',num2str(IDXlist(4))))


%%%%%%%%%%%%%%%%%%%%%%%%
%以下、local関数(getinput, cal_x2x)
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%getinput:実験ログ（T）からshot番号（IDX）におけるオペレーションの値を出力
%出力は構造体に変更しても良いかも
function [date, shot, TF_shot, offset_TF, i_EF, start, v_TF, d_tacq, d_tacqTF, trange, n] = getinput(T,IDX)
date=T.date(IDX);
shot=T.shot(IDX);
TF_shot=T.TFoffset(IDX);
offset_TF=isfinite(TF_shot);

if isnan(T.EF_A_(IDX))%%NaNでないことを確認（ログが空白だとNaNになる）
    i_EF=150;
else  %NaNなら150をとりあえず代入、記入されているときはその値を使う
    i_EF=T.EF_A_(IDX);
end

start=T.Period_StartTime_(IDX);
v_TF=T.TF_kV_(IDX); %TF[kV]

d_tacq=T.d_tacq(IDX);
d_tacqTF=T.TFdtacq(IDX);

trange=470:500;

n=50; %rz方向のメッシュ数
end

%%%plot_psi:pcbプローブの磁気面・電流密度・X点・O点の時系列プロット
function [xEt,xJt,eta,trange]=def_sheet(T, IDX,pathname)
[date, ~, ~, ~, i_EF, ~, ~, d_tacq, d_tacqTF, trange, n] = getinput(T,IDX);%Tのテーブルから入力のリストを出力
tsize = numel(trange);
cal_trange=trange(1):trange(tsize)+1;

[grid2D, data2D] = pcbdata(date, d_tacq, d_tacqTF, cal_trange, [], n,i_EF);

if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end

% %%%速度Vの各成分を計算
% Vz=zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2));
% Vr=zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2));
% for i=1:size(trange,2)-1
%     Vz(:,:,i)=data2D.Et(:,:,i)./data2D.Br(:,:,i);
%     Vr(:,:,i)=-data2D.Et(:,:,i)./data2D.Bz(:,:,i);
% end

%%%midplaneとかO点、X点を探す
[psimid,mid]=min(data2D.psi,[],2);
[opoint,~]=islocalmin(psimid,1);
[xpoint,~]=islocalmax(psimid,1);
[xp_psi,maxxp]=max(squeeze(psimid),[],1);
% onum=squeeze(sum(opoint,1));
r_index=1:n;

%%%X点での実効抵抗値、Et,Jt
% Et_midmax=zeros(size(trange,2)-1,1);
% Jt_midmin=zeros(size(trange,2)-1,1);
xEt=zeros(tsize,1);
xJt=zeros(tsize,1);
eta=zeros(tsize,1);
% Et_mid=zeros(n,1,tsize);
% Jt_mid=zeros(n,1,tsize);
% mid_logical=false(n,n,tsize);

for i=1:tsize
    if sum(maxxp(:,i))>0
        r_maxxp=r_index(maxxp(:,i));
        r_s=max(1,r_maxxp-2);
        r_n=min(n,r_maxxp+2);
        z_maxxp=squeeze(mid(maxxp(:,i),:,i));
        z_w=max(1,z_maxxp-2);
        z_e=min(n,z_maxxp+2);
        xEt(i,1)=sum(data2D.Et(r_s:r_n,z_w:z_e,i),'all')/((r_n-r_s+1)*(z_e-z_w+1));
        xJt(i,1)=sum(data2D.Jt(r_s:r_n,z_w:z_e,i),'all')/((r_n-r_s+1)*(z_e-z_w+1));
        eta(i,1)=xEt(i,1)./xJt(i,1);
    else
        xEt(i,1)=NaN;
        xJt(i,1)=NaN;
        eta(i,1)=NaN;
    end
end

Etind=isoutlier(xEt,'median',1);
xEt(Etind)=NaN;
Jtind=isoutlier(xJt,'median',1);
xJt(Jtind)=NaN;
etaind=isoutlier(eta,'median',1);
eta(etaind)=NaN;


%X点でのJtのプロット
figure('Position', [0 0 700 700],'visible','on')
subplot(3,1,1)
plot(trange, xEt,'--o')
xlabel('Time [us]')
ylabel('Et at X[V/m]')
xlim([trange(1) trange(tsize)])

subplot(3,1,2)
plot(trange, xJt,'--o')
xlabel('Time [us]')
ylabel('Jt at X[A/m^{2}]')
xlim([trange(1) trange(tsize)])

subplot(3,1,3)
plot(trange, eta,'--o')
xlabel('Time [us]')
ylabel('Resistivity [Ω]')
xlim([trange(1) trange(tsize)])

sgtitle(strcat('IDX=',num2str(IDX)))

% saveas(gcf,strcat(pathname.save,'\x_IDX',num2str(IDX),'.png'))
% close

% for i=1:tsize
%     for j=1:n
%         mid_logical(j,mid(j,:,i),i)=true;
%     end
%     Et_mid(:,:,i)=data2D.Et(mid_logical(:,:,i));
%     Jt_mid(:,:,i)=data2D.Jt(mid_logical(:,:,i));
% %     Et_midmax(i)=max(Et_mid(:,:,i),[],'all');
% %     Jt_midmin(i)=min(Jt_mid(:,:,i),[],'all');
%     if sum(xpoint(:,1,i))>0
%         xEt(i)=Et_mid(maxxp(:,i),:,i);
%         xJt(i)=Jt_mid(maxxp(:,i),:,i);
%     end
% %     eta_all=0;
% %     if sum(xpoint(:,1,i))>0
% %         eta_all=Et_mid(xpoint(:,:,i),:,i)./Jt_mid(xpoint(:,:,i),:,i);
% %     end
% %     eta(i)=mean(eta_all); %X点が複数個ある場合は平均値を採用（maxでも良いかも）
% end

% figure
% plot(trange,eta(:,1))
% title(strcat('IDX=',num2str(IDX)))
% xlabel('time [us]')
% ylabel('Resistivity')

%%%current sheet
%sheet=data2D.Jt



% f=figure;
% f.WindowState = 'maximized';
%  start=11; %460+?
%  t_start=460+start;
%  for m=1:10 %図示する時間
%      i=start+m; %end
%      t=trange(i);
%      subplot(2,5,m)
%     contourf(grid2D.zq(1,:),grid2D.rq(:,1),v_z(:,:,i),20,'LineStyle','none')
%     colormap(jet)
%     axis image
%     axis tight manual
%     %     xlim([-0.02 0.02])
%     %     ylim([0.12 0.27])
%     %caxis([-2.7*1e+6,3*1e+6]) %カラーバーの軸の範囲
%     %caxis([-maxrange,maxrange])
%     colorbar('Location','eastoutside')
%     %zlim([-1 1])
%     %colormap(bone)
%     %%カラーバーのラベル付け
%     %c = colorbar;
%     %c.Label.String = 'Jt';
%     hold on
%     %plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
%     %contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),35,'black')
%     plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"ro")
%     plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"rx")
%     hold off
%     title(string(t)+'us')
%     xlabel('z')
%     ylabel('r')
%  end
 
 %sgtitle(strcat('IDX=',num2str(IDX),': shot=',num2str(date),num2str(shot,'%03i'),': dtacq=',num2str(T.d_tacq(IDX))))

end

