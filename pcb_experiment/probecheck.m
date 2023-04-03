%%%%%%%%%%%%%%%%%%%%%%%%
%pcbプローブと装置の磁場信号極性チェック
%dtacqのshot番号を直接指定する場合
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%ここが各PCのパス
%【※コードを使用する前に】環境変数を設定しておくか、matlab内のコマンドからsetenv('パス名','アドレス')で指定してから動かす
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.save=getenv('savedata_path');%outputデータ保存先

pathname.rawdata=getenv('rawdata_path');%dtacqのrawdataの保管場所
pathname.rawdata38=getenv('rawdata038_path');%dtacq a038のrawdataの保管場所
pathname.woTFdata=getenv('woTFdata_path');%rawdata（TFoffset引いた）の保管場所

%%%%実験オペレーションの取得
%直接入力の場合
dtacqlist=10314; %【input】dtacqの保存番号
date = 230103;%【input】計測日
d_tacqTF = 10312;%【input】TFoffsetのdtacq保存番号

%磁気面出す場合は適切な値を入力、磁場信号のみプロットする場合は変更不要
i_EF = 0;%【input】EF電流
trange=460:490;%【input】計算時間範囲
n=50; %【input】rz方向のメッシュ数

% spread sheetから ログのテーブルを取得してTに格納、参照する場合
% %Github/test-open/getTS6log.mを使用
% DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
% T=getTS6log(DOCID);
%  [date, shot, TF_shot, offset_TF, i_EF, start, Doppler_t, d_tacq, d_tacqTF, trange, t, n] = getinput(T,IDX)

for d_tacq=dtacqlist(1,1)
check_signal(date, d_tacq, d_tacqTF,trange, n, i_EF, pathname); %通常の時系列プロット
end

%%%%%%%%%%%%%%%%%%%%%%%%
%以下、local関数(getinput, plot_psi)
%%%%%%%%%%%%%%%%%%%%%%%%

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

trange=300:600;
t=T.DopplerDelay(IDX);
n=50; %rz方向のメッシュ数
end


function check_signal(date, d_tacq, d_tacqTF,trange, n, i_EF, pathname)
% filename=strcat(pathname.rawdata38,'rawdata_noTF_dtacq',num2str(d_tacq),'.mat');
% load(filename,'rawdata');
filename=strcat(pathname.rawdata,'\rawdata_dtacq38_shot',num2str(d_tacq),'_tfshot',num2str(d_tacqTF),'.mat');
load(filename,'rawdata');%1000×192

%正しくデータ取得できていない場合はreturn
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
y_upper_lim = 0.05;%3e-3;%0.1;%縦軸プロット領域（b_z上限）
y_lower_lim = -0.05;%3e-3;%-0.1;%縦軸プロット領域（b_z下限）
t_start=350;%430;%455;%横軸プロット領域（開始時間）
t_end=600;%550;%横軸プロット領域（終了時間）
r_ch=col1+col2;%r方向から挿入した各プローブのチャンネル数

f=figure;
f.WindowState = 'maximized';
for i=1:r
    for j=1:col1
        subplot(r,col1,(i-1)*col1+j)
        if ok(r_ch*(i-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bz(t_start:t_end,r_ch*(i-1)+j))
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
% saveas(gcf,strcat(pathname.save,'\date',num2str(date),'_dtacq',num2str(d_tacq),'_01','.png'))
% close

f2=figure;
f2.WindowState = 'maximized';
for i=1:r
    for j=col1+1:col1+col2
        subplot(r,col2,(i-1)*col2+j-col1)
        if ok(r_ch*(i-1)+j)==1 
            plot(t_start:t_end,bz(t_start:t_end,r_ch*(i-1)+j))
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
% saveas(gcf,strcat(pathname.save,'\date',num2str(date),'_dtacq',num2str(d_tacq),'_02','.png'))
% close


end
