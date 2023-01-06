%%%%%%%%%%%%%%%%%%%%%%%%
%200ch用新規pcbプローブと装置の磁場信号極性チェック
%dtacqのshot番号を直接指定する場合
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%ここが各PCのパス
%【※コードを使用する前に】環境変数を設定しておくか、matlab内のコマンドからsetenv('パス名','アドレス')で指定してから動かす
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.save=getenv('savedata_path');%outputデータ保存先

pathname.rawdata38=getenv('rawdata038_path');%dtacq a038のrawdataの保管場所
pathname.woTFdata=getenv('woTFdata_path');%rawdata（TFoffset引いた）の保管場所

pathname.rawdata=getenv('rawdata_path');%dtacqのrawdataの保管場所

%%%%実験オペレーションの取得
%直接入力の場合
dtacqlist=39;
shotlist=240;%241;%【input】dtacqの保存番号
tfshotlist=0;%237;
date = 221223;%【input】計測日
n=numel(shotlist);%計測データ数

% %磁気面出す場合は適切な値を入力、磁場信号のみプロットする場合は変更不要
% d_tacqTF = '';%【input】TFoffsetのdtacq保存番号
% i_EF = 0;%【input】EF電流
% trange=460:490;%【input】計算時間範囲
% n=50; %【input】rz方向のメッシュ数

for i=1:n
    dtacq_num=dtacqlist(i);
    shot=shotlist(i);
    tfshot=tfshotlist(i);
    check_signal(date, dtacq_num, shot, tfshot, pathname); 
end

%%%%%%%%%%%%%%%%%%%%%%%%
%以下、local関数
%%%%%%%%%%%%%%%%%%%%%%%%

function check_signal(date, dtacq_num, shot, tfshot, pathname)
filename=strcat(pathname.rawdata,'\rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat');
% filename=strcat(pathname.rawdata,'rawdata_noTF_dtacq',num2str(d_tacq),'.mat');
load(filename,'rawdata');%1000×192

%正しくデータ取得できていない場合はreturn
if numel(rawdata)< 500
    return
end

%較正係数のバージョンを日付で判別
sheets = sheetnames('coeff200ch.xlsx');
sheets = str2double(sheets);
sheet_date=max(sheets(sheets<=date));

C = readmatrix('coeff200ch.xlsx','Sheet',num2str(sheet_date));
ok = logical(C(:,14));
P=C(:,13);
coeff=C(:,12);
zpos=C(:,9);
rpos=C(:,10);
probe_num=C(:,5);
probe_ch=C(:,6);
ch=C(:,7);
d2p=C(:,15);
d2bz=C(:,16);
d2bt=C(:,17);

b=rawdata.*coeff';%較正係数RC/NS
b=b.*P';%極性揃え
b=smoothdata(b,1);

%デジタイザchからプローブ通し番号順への変換
bz=zeros(1000,100);
bt=bz;
ok_bz=true(1,100);
ok_bt=ok_bz;

for i=1:192
    if rem(ch(i),2)==1
        bz(:,ceil(ch(i)/2))=b(:,i);
        ok_bz(ceil(ch(i)/2))=ok(i);
    elseif rem(ch(i),2)==0
        bt(:,ch(i)/2)=b(:,i);
        ok_bt(ceil(ch(i)/2))=ok(i);
    end
end

Pcheck=[1	-1	1	1	-1	-1	-1	1	1	1	1	1	1	-1	-1	0	-1	-1	1	-1	1	0	1	-1	-1	1	-1	-1	1	-1	-1	1	-1	1	-1	1	-1	-1	-1	1	-1	-1	1	1	-1	-1	-1	-1	-1	-1	1	-1	-1	-1	-1	-1	-1	-1	-1	-1	1	1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	1	-1	1	-1	-1	-1	-1	1	-1	-1	1	-1	-1	1	-1	-1	1	1	1	-1	1	-1	-1	1	1	1	-1];
bz=bz.*Pcheck;
ok_bz([5 6 11 16 20 22 31 33 39 49 63 66 71 72 79 80 95 100])=false;

% log_d2bz=~isnan(d2bz);
% ok_bz=ok(log_d2bz);
% d2bz=d2bz(log_d2bz);
% bz=b(:,d2bz');
% ok_bz=ok_bz(d2bz);
% 
% log_d2bt=~isnan(d2bt);
% ok_bt=ok(log_d2bt);
% d2bt=d2bt(log_d2bt);
% bt=b(:,d2bt');
% ok_bt=ok_bt(d2bt);
% 
% log_d2p=~isnan(d2p);
% ok=ok(log_d2p);
% d2p=d2p(log_d2p);
% b=b(:,d2p');
% ok=ok(d2p);

%生信号描画用パラメータ
r = 5;%プローブ本数＝グラフ出力時の縦に並べる個数
col = 10;%グラフ出力時の横に並べる個数
y_upper_lim = 0.02;%3e-3;%0.1;%縦軸プロット領域（b_z上限）
y_lower_lim = -0.02;%3e-3;%-0.1;%縦軸プロット領域（b_z下限）
t_start=300;%430;%455;%横軸プロット領域（開始時間）
t_end=600;%550;%横軸プロット領域（終了時間）
% r_ch=col1+col2;%r方向から挿入した各プローブのチャンネル数

f1=figure;
f1.WindowState = 'maximized';
for i=1:r
    for j=1:col
        subplot(r,col,(i-1)*col+j)
        if ok_bz(col*(i-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bz(t_start:t_end,col*(i-1)+j))
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bz(t_start:t_end,col*(i-1)+j),'r:')
        end   
        title(num2str(2.*(col*(i-1)+j)-1));
        xticks([t_start t_end]);
        ylim([y_lower_lim y_upper_lim]);
    end
end
sgtitle('Bz signal probe1-5')

f2=figure;
f2.WindowState = 'maximized';
for i=1:r
    for j=1:col
        subplot(r,col,(i-1)*col+j)
        if ok_bz(col*(i+r-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bz(t_start:t_end,col*(i+r-1)+j))
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bz(t_start:t_end,col*(i+r-1)+j),'r:')
        end   
        title(num2str(2.*(col*(i+r-1)+j)-1));
        xticks([t_start t_end]);
        ylim([y_lower_lim y_upper_lim]);
    end
end
sgtitle('Bz signal probe6-10')

% f3=figure;
% f3.WindowState = 'maximized';
% for i=1:r
%     for j=1:col
%         subplot(r,col,(i-1)*col+j)
%         if ok_bt(col*(i-1)+j)==1 %okなチャンネルはそのままプロット
%             plot(t_start:t_end,bt(t_start:t_end,col*(i-1)+j))
%         else %NGなチャンネルは赤色点線でプロット
%             plot(t_start:t_end,bt(t_start:t_end,col*(i-1)+j),'r:')
%         end   
%         title(num2str(2.*(col*(i-1)+j)));
%         xticks([t_start t_end]);
%         ylim([y_lower_lim y_upper_lim]);
%     end
% end
% sgtitle('Bt signal probe1-5')
% 
% f4=figure;
% f4.WindowState = 'maximized';
% for i=1:r
%     for j=1:col
%         subplot(r,col,(i-1)*col+j)
%         if ok_bt(col*(i+r-1)+j)==1 %okなチャンネルはそのままプロット
%             plot(t_start:t_end,bt(t_start:t_end,col*(i+r-1)+j))
%         else %NGなチャンネルは赤色点線でプロット
%             plot(t_start:t_end,bt(t_start:t_end,col*(i+r-1)+j),'r:')
%         end   
%         title(num2str(2.*(col*(i+r-1)+j)));
%         xticks([t_start t_end]);
%         ylim([y_lower_lim y_upper_lim]);
%     end
% end
% sgtitle('Bt signal probe6-10')

% saveas(gcf,strcat(pathname.save,'\date',num2str(date),'_dtacq',num2str(d_tacq),'_02','.png'))
% close

end
