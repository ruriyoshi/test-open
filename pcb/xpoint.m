%%%%%%%%%%%%%%%%%%%%%%%%
%　Xポイントの座標
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
pathname.save='C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\過去データ解析\220928'; %保存先
%pathname.rawdata='/Users/mgar/rawdata_a038/'; %rawdataの保管場所

%%%%(2)ログから解析したいデータを検索
%Github/test-open/searchlog.mを使用

% node='date';  % 【input】検索する列の名前. T.Properties.VariableNamesで一覧表示できる
%  pat=211223;   % 【input】検索パターン（数値なら一致検索、文字なら含む検索）　
% 
% searchlog(T,node,pat); % ログのテーブルから当てはまるものを抽出した新しいテーブルを作成

%%%%(3)指定したshotの解析
% IDXlist=[2897 2906 2907 2912 2913] ; %2870:2921; %【input】テーブルから解析したいshot番号を抽出して入力
IDXlist=[2911:2913 2925 2926 2927 2931 2933 2947:2950 2942 2943 2946];
for IDX=IDXlist(1,1)
[xr,xz]=plot_xpoint(T, pathname,IDX);
end

%%%%%%%%%%%%%%%%%%%%%%%%
%以下、local関数
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%getinput:実験ログ（T）からshot番号（IDX）におけるオペレーションの値を出力
%出力は構造体に変更しても良いかも
function [date, shot, TF_shot, offset_TF, i_EF, start, Doppler_t, d_tacq, d_tacqTF, trange, t, n] = getinput(T,IDX)
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
Doppler_t=T.DopplerDelay(IDX);

d_tacq=T.d_tacq(IDX);
d_tacqTF=T.TFdtacq(IDX);

trange=460:490;
t=T.DopplerDelay(IDX);
n=50; %rz方向のメッシュ数
end

function [xr,xz]=plot_xpoint(T, pathname,IDX);
[date, shot, TF_shot, offset_TF, i_EF, start, Doppler_t, d_tacq, d_tacqTF,trange, t, n] = getinput(T,IDX);%Tのテーブルから入力のリストを出力

[grid2D, data2D] = pcbdata(date, d_tacq,d_tacqTF,trange, [], n,i_EF);

if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end
    maxrange=max(abs(data2D.Jt),[],'all');

%%%midplaneとかO点、X点を探す
[psimid,mid]=min(data2D.psi,[],2);
[opoint,p]=islocalmin(psimid,1);
[xpoint,~]=islocalmax(psimid,1);
[xp_psi,maxxp]=max(squeeze(psimid),[],1);

tsize=numel(trange);
xr=zeros(1,tsize);%X-point(1点)のr座標index
xz=xr;%X-point(1点)のz座標index
r_index=1:n;

for i=1:tsize
    %(1)search Xpoint posision (index形式・1点のみ)：(z_xind,r_xind)
    if sum(xpoint(:,:,i))>0
        r_xind0=r_index(xpoint(:,:,i));
        z_xind0=mid(xpoint(:,:,i),:,i);
        if numel(r_xind0)>1
            [~,I]=max(psimid(r_xind0));
            r_xind0=r_xind0(I);
            z_xind0=z_xind0(I);
        end
        xr(1,i)=grid2D.rq(r_xind0,1);
        xz(1,i)=grid2D.zq(1,z_xind0);
        clear z_xind0 r_xind0
        
    else
        xr(1,i)=NaN;
        xz(1,i)=NaN;
    end
end

%X点r座標プロット
figure
plot(trange,xr,'b*','MarkerSize',8)
ylabel('r [m]')
xlabel('Time [us]')
% xlim([460 482])
title(strcat('IDX=',num2str(IDX)))
ha1 = gca;
ha1.LineWidth = 1;
ha1.FontSize=13;
% saveas(gcf,strcat(pathname.save,'\',num2str(IDX),'_xr.png'))
% close

% % %%磁気面時間発展プロット
% f=figure;
% f.WindowState = 'maximized';
%  start=20; %460+?
%  for m=1:10 %図示する時間
%      i=start+m; %end
%      t=trange(i);
%      subplot(2,5,m)
%     contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),10,'LineStyle','none')
%     colormap(jet)
%     axis image
%     axis tight manual
% %     caxis([-2.5*1e+6,2.5*1e+6]) %カラーバーの軸の範囲
%     %caxis([-1*1e+6,1*1e+6])
%     %caxis([-maxrange,maxrange])
%     colorbar('Location','eastoutside')
%     %カラーバーのラベル付け
% %     c = colorbar;
% %     c.Label.String = 'Jt [A/m^{2}]';
%     hold on
%     plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
%     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),50,'black')
%     plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"bo")
%     plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"bx")
%     hold off
%     title(string(t)+' us')
%     xlabel('z [m]')
%     ylabel('r [m]')
%  end
%  sgtitle(strcat('IDX=',num2str(IDX)))

end
