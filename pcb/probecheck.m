%%%%%%%%%%%%%%%%%%%%%%%%
%pcbプローブと装置の極性チェック
%dtacqのshot番号を直接指定する場合
%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%(1)spread sheetから ログのテーブルを取得してTに格納
% %Github/test-open/getTS6log.mを使用
% DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
% T=getTS6log(DOCID);

%%%%%ここが各PCのパス
%環境変数を設定していない場合はパスを''内に全て記入する（使用しないパスは空白''で良い）
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.save='C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\221017'; %保存先

pathname.rawdata='C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\rawdata_a038_noTF\'; %rawdataの保管場所

%%%%(3)指定したshotの解析
dtacqlist=3327;%10002:10004; %【input】テーブルから解析したいshot番号を抽出して入力

date = 211214;%221014;
d_tacqTF = '';
i_EF = 0;
trange=460:490;
n=50; %rz方向のメッシュ数

for d_tacq=dtacqlist(1,1)
plot_psi(date, d_tacq, d_tacqTF,trange, n, i_EF, pathname); %通常の時系列プロット
%plot_position(date, d_tacq, d_tacqTF,trange, n, i_EF, pathname); %計測位置、各位置での生信号も含めた確認用プロット
end

%%%%%%%%%%%%%%%%%%%%%%%%
%以下、local関数(getinput, plot_psi, plot_posision)
%%%%%%%%%%%%%%%%%%%%%%%%

%%%plot_psi:pcbプローブの磁気面・電流密度・X点・O点の時系列プロット
%Github/test-open/pcbdata.mを使用->[チャンネルごとの生信号のプロット]の項目をコメントアウトしなければ、別figureで確認用にプロットされる
function plot_psi(date, d_tacq, d_tacqTF,trange, n, i_EF, pathname)

filename=strcat(pathname.rawdata,'rawdata_noTF_dtacq',num2str(d_tacq),'.mat');
load(filename,'rawdata');

if numel(rawdata)< 500
    grid2D=NaN;
    data2D=NaN;
    return
end

load('rc_coeff2020.mat')
    
[ok, bz, rpos, zpos,p_ch] = getpcbbz(rawdata, coeff,date);

% チャンネルごとの生信号のプロット
bz=smoothdata(bz,1);

%生信号描画用パラメータ
r = 5;%プローブ本数＝グラフ出力時の縦に並べる個数
col1 = 12;%1枚目のグラフ出力時の横に並べる個数
col2 = 13;%2枚目のグラフ出力時の横に並べる個数
y_upper_lim = inf;%0.1;%縦軸プロット領域（b_z上限）
y_lower_lim = -inf;%-0.1;%縦軸プロット領域（b_z下限）
t_start=430;%455;%横軸プロット領域（開始時間）
t_end=520;%横軸プロット領域（終了時間）
r_ch=col1+col2;%r方向から挿入した各プローブのチャンネル数


f=figure;
f.WindowState = 'maximized';
for i=1:r
    for j=1:col1
        subplot(r,col1,(i-1)*col1+j)
        if ok(r_ch*(i-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bz(t_start:t_end,r_ch*(i-1)+j))
            %plot(t_start:t_end,bz_s(t_start:t_end,r_ch*(i-1)+j),'r')
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bz(t_start:t_end,r_ch*(i-1)+j),'r:')
        end   
        title(num2str(p_ch(i,j)));
        %xlim([t_start t_end]);
        xticks([t_start t_end]);
        ylim([y_lower_lim y_upper_lim]);
        %ylim([-0.02 0.04]);
    end
end
saveas(gcf,strcat(pathname.save,'\date',num2str(date),'_dtacq',num2str(d_tacq),'_01','.png'))

f2=figure;
f2.WindowState = 'maximized';
for i=1:r
    for j=col1+1:col1+col2
        subplot(r,col2,(i-1)*col2+j-col1)
        if ok(r_ch*(i-1)+j)==1 
            plot(t_start:t_end,bz(t_start:t_end,r_ch*(i-1)+j))
            %plot(t_start:t_end,bz_s(t_start:t_end,r_ch*(i-1)+j),'r')
        else 
            plot(t_start:t_end,bz(t_start:t_end,r_ch*(i-1)+j),'r:')
        end   
        title(num2str(p_ch(i,j)));
        %xlim([t_start t_end]);
        xticks([t_start t_end]);
        ylim([y_lower_lim y_upper_lim]);
        %ylim([-0.02 0.04]);
    end
end
saveas(gcf,strcat(pathname.save,'\date',num2str(date),'_dtacq',num2str(d_tacq),'_02','.png'))

% %磁気面描画
% [zq,rq]=meshgrid(linspace(min(zpos),max(zpos),n),linspace(min(rpos),max(rpos),n));
% grid2D=struct('zq',zq,'rq',rq);
% clear zq rq
% 
% data2D = data2Dcalc(EF, grid2D, n, trange, rpos, zpos, bz, ok);
% 
% [grid2D, data2D] = pcbdata(date, d_tacq, d_tacqTF,trange, [], n, i_EF);
% 
% if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
%     return
% end
%     maxrange=max(abs(data2D.Jt),[],'all');
% 
% %%%midplaneとかO点、X点を探す
% [psimid,mid]=min(data2D.psi,[],2); %各r,tでのpsiの最小値,時間
% [opoint,p]=islocalmin(psimid,1); %全rでのpsiの極小値
% [xpoint,~]=islocalmax(psimid,1); %全rでのpsiの極大値
% [xp_psi,maxxp]=max(squeeze(psimid),[],1);
% % onum=squeeze(sum(opoint,1));
% % trange(onum~=0)
%     %maxrange=2e6;

% %%磁気面時間発展プロット
% f=figure;
% f.WindowState = 'maximized';
%  start=0; %460+?
%  t_start=460+start;
%  for m=1:10 %図示する時間
%      i=start+m; %end
%      t=trange(i);
%      subplot(2,5,m)
%     contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Jt(:,:,i),10,'LineStyle','none')
%     colormap(jet) %jet/parula
%     axis image
%     axis tight manual
%     %     xlim([-0.02 0.02])
%     %     ylim([0.12 0.27])
%     caxis([-10*1e+6,10*1e+6]) %カラーバーの軸の範囲
%     %caxis([-maxrange,maxrange])
%     colorbar('Location','eastoutside')
%     %zlim([-1 1])
%     %colormap(bone)
%     %%カラーバーのラベル付け
%     %c = colorbar;
%     %c.Label.String = 'Jt';
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
% 
% sgtitle(strcat('date=',num2str(date),': dtacq=',num2str(d_tacq)));
% 
%figureの保存
%saveas(gcf,strcat(pathname.save,'\date',num2str(date),'_dtacq',num2str(d_tacq),'_time',num2str(t_start),'.png'))
%save(strcat(pathname.save,'\date',num2str(date),'_dtacq',num2str(d_tacq),'.mat'))

%close
end
