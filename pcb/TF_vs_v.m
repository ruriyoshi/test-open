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
% pathname.ts3u='ts3u_path';%old-koalaのts-3uまでのパス（mrdなど）
% pathname.fourier='fourier_path';%fourierのmd0（データックのショットが入ってる）までのpath
% pathname.NIFS='NIFS_path';%resultsまでのpath（ドップラー、SXR）
% pathname.save='/Users/mgar/pcb_save'; %保存先

%%%%(2)ログから解析したいデータを検索
%Github/test-open/searchlog.mを使用

% node='date';  % 【input】検索する列の名前. T.Properties.VariableNamesで一覧表示できる
%  pat=211223;   % 【input】検索パターン（数値なら一致検索、文字なら含む検索）　
% 
% searchlog(T,node,pat); % ログのテーブルから当てはまるものを抽出した新しいテーブルを作成

%%%%(3)指定したshotの解析
IDXlist=2870:2920; %【input】テーブルから解析したいshot番号を抽出して入力
for IDX=IDXlist(1,37)
cal_v(T,IDX); 
end

% figure
% for i=1:size(IDXlist,2)
%     IDX=IDXlist(1,i);
%     [x2x,v_TF,~]=cal_x2x(T,IDX);
%     if x2x>0
%         plot(v_TF,x2x,"ro");
%     end
%     clear v_TF x2x
%     hold on
% end
% hold off
% title(strcat('IDX=',num2str(IDXlist(1,1)),':',num2str(IDXlist(1,i))))
% xlabel('TF [kV]')
% ylabel('X to X [m]')
% 
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

trange=460:490;

n=50; %rz方向のメッシュ数
end

%%%plot_psi:pcbプローブの磁気面・電流密度・X点・O点の時系列プロット
function cal_v(T,IDX)
[date, ~, ~, ~, i_EF, ~, v_TF, d_tacq, d_tacqTF, trange, n] = getinput(T,IDX);%Tのテーブルから入力のリストを出力

[grid2D, data2D] = pcbdata(date, d_tacq, d_tacqTF, trange, [], n,i_EF);

if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end

%速度Vの各成分を計算
Vz=zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2));
Vr=zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2));
for i=1:size(trange,2)-1
    Vz(:,:,i)=data2D.Et(:,:,i)./data2D.Br(:,:,i);
    Vr(:,:,i)=-data2D.Et(:,:,i)./data2D.Bz(:,:,i);
end


%%%midplaneとかO点、X点を探す
[psimid,mid]=min(data2D.psi,[],2);
[opoint,~]=islocalmin(psimid,1);
[xpoint,~]=islocalmax(psimid,1);
%[xp_psi,maxxp]=max(squeeze(psimid),[],1);
% onum=squeeze(sum(opoint,1));

%%%X点での実効抵抗値
eta=zeros(size(trange,2),1);
Et_mid=zeros(n,1,size(trange,2)-1);
Jt_mid=zeros(n,1,size(trange,2)-1);
mid_logical=false(n,n,size(trange,2)-1);
for i=1:size(trange,2)-1
    for j=1:n
        mid_logical(j,mid(j,:,i),i)=true;
    end
    Et_mid(:,:,i)=data2D.Et(mid_logical(:,:,i));
    Jt_mid(:,:,i)=data2D.Jt(mid_logical(:,:,i));
    eta_all=0;
    if sum(xpoint(:,1,i))>0
        eta_all=Et_mid(xpoint(:,:,i),:,i)./Jt_mid(xpoint(:,:,i),:,i);
    end
    eta(i)=mean(eta_all); %X点が複数個ある場合は平均値を採用（maxでも良いかも）
end

figure
plot(trange,eta(:,1))
title(strcat('IDX=',num2str(IDX)))
xlabel('time [us]')
ylabel('Resistivity')



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

