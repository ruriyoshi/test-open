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
dtacqlist=40;
shotlist=13;%【input】dtacqの保存番号
tfshotlist=0;
date = 230520;%【input】計測日
n=numel(shotlist);%計測データ数

trange = 400:600;

% %磁気面出す場合は適切な値を入力、磁場信号のみプロットする場合は変更不要
% d_tacqTF = '';%【input】TFoffsetのdtacq保存番号
% i_EF = 0;%【input】EF電流
% trange=460:490;%【input】計算時間範囲
% n=50; %【input】rz方向のメッシュ数

for i=1:n
    dtacq_num=dtacqlist(i);
    shot=shotlist(i);
    tfshot=tfshotlist(i);
    low_n_mode(date, dtacq_num, shot, tfshot, pathname,trange); 
end

%%%%%%%%%%%%%%%%%%%%%%%%
%以下、local関数
%%%%%%%%%%%%%%%%%%%%%%%%

function low_n_mode(date, dtacq_num, shot, tfshot, pathname,trange)
filename=strcat(pathname.rawdata,'rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat');
% filename=strcat(pathname.rawdata,'rawdata_noTF_dtacq',num2str(d_tacq),'.mat');
load(filename,'rawdata');%1000×192

%正しくデータ取得できていない場合はreturn
if numel(rawdata)< 500
    return
end

%較正係数のバージョンを日付で判別
sheets = sheetnames('C:\Users\uswk0\OneDrive\デスクトップ\Github\test-open\pcb_experiment_new\low_n_a040.xlsx');
sheets = str2double(sheets);
sheet_date=max(sheets(sheets<=date));

C = readmatrix('C:\Users\uswk0\OneDrive\デスクトップ\Github\test-open\pcb_experiment_new\low_n_a040.xlsx','Sheet',num2str(sheet_date));
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

b=rawdata(:,17:40).*coeff';%較正係数RC/NS
b=b.*P';%極性揃え
b=smoothdata(b,1);

%デジタイザchからプローブ通し番号順への変換
bz=zeros(1000,8);
bt=bz;
br=bz;
ok_bz=true(1,8);
ok_bt=ok_bz;

for i=1:24
    if probe_ch(i)==1
        bz(:,probe_num(i))=b(:,i);
        
    elseif probe_ch(i)==2
        bt(:,probe_num(i))=b(:,i);

    elseif probe_ch(i)==3
        br(:,probe_num(i))=b(:,i);
       
    end
end


    
x = [0,4/24,7/24,9/24,12/24,15/24,17/24,20/24]; % non-uniform probe placement

t = trange;

%measured_dBzdt = measured_dBzdt(3500:6000,:);
[n_Amp,n_Ph] = toroidal_mode(t,x,bz); % get amplitude and phase
% calculate frequency from phase (smoothed)
n1_Omega = calculate_omega(smooth(cumulative_phase(n_Ph(1,:)),0.1,'loess'));
n2_Omega = calculate_omega(smooth(cumulative_phase(n_Ph(2,:)),0.1,'loess'));
n3_Omega = calculate_omega(smooth(cumulative_phase(n_Ph(3,:)),0.1,'loess'));


%*************** plotting ***************
figure('Position', [0 0 800 800]);
subplot(2,1,1)
hold on
plot(t,n_Amp(1,:));
plot(t,n_Amp(2,:),'k');
plot(t,n_Amp(3,:),'b');
plot(t,n_Amp(4,:),'r');
legend('n=0','n=1','n=2','n=3');
ylim([0 2000])
xlim([400 600])
ylabel('[T/s]')
xlabel('[us]')
hold off

subplot(2,1,2)
hold on
%plot(t,smooth(n_Amp(2,:)./n_Amp(1,:),0.08,'loess'),'k');
%plot(t,smooth(n_Amp(3,:)./n_Amp(1,:),0.08,'loess'),'b');
%plot(t,smooth(n_Amp(4,:)./n_Amp(1,:),0.08,'loess'),'r');

plot(t,n_Amp(2,:)./n_Amp(1,:),'k');
plot(t,n_Amp(3,:)./n_Amp(1,:),'b');
plot(t,n_Amp(4,:)./n_Amp(1,:),'r');
legend('n=1/n=0','n=2/n=0','n=3/n=0');

%legend('n=1','n=2','n=3');
%ylim([0 1]);
%ylim([0 8000])
xlim([400 600])
xlabel('[us]')

hold off
end
