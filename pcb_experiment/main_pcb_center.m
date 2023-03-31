%%%%%%%%%%%%%%%%%%%%%%%%
%125ch秋光プローブ＋200ch用新規pcbプローブのみでの磁気面（Bz）
%dtacqのshot番号を直接指定する場合
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%ここが各PCのパス
%【※コードを使用する前に】環境変数を設定しておくか、matlab内のコマンドからsetenv('パス名','アドレス')で指定してから動かす
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.save=getenv('savedata_path');%outputデータ保存先

pathname.rawdata038=getenv('rawdata038_path');%dtacq a038のrawdataの保管場所
pathname.woTFdata=getenv('woTFdata_path');%rawdata（TFoffset引いた）の保管場所

pathname.rawdata=getenv('rawdata_path');%dtacqのrawdataの保管場所

%%%%実験オペレーションの取得
%直接入力の場合
dtacqlist=[38 39];
shotlist=[10694 632];%【input】dtacqの保存番号
tfshotlist=[0 0];
date = 230127;%【input】計測日
n_data=numel(shotlist(:,1));%計測データ数

i_EF = 0;%【input】EF電流
trange=450:510;%【input】計算時間範囲
n=51; %【input】rz方向のメッシュ数

for i=1:n_data
    dtacq_num=dtacqlist(i,:);
    shot=shotlist(i,:);
    tfshot=tfshotlist(i,:);
    plot_psi_center(date, dtacq_num, shot, tfshot, pathname,n,i_EF,trange); 
end

%%%%%%%%%%%%%%%%%%%%%%%%
%以下、local関数
%%%%%%%%%%%%%%%%%%%%%%%%

function plot_psi_center(date, dtacq_num, shot, tfshot, pathname, n,i_EF,trange)

filename1=strcat(pathname.rawdata,'rawdata_dtacq',num2str(dtacq_num(1)),'_shot',num2str(shot(1)),'_tfshot',num2str(tfshot(1)),'.mat');
load(filename1,'rawdata');%a038
rawdata1=rawdata;
clear rawdata
filename2=strcat(pathname.rawdata,'rawdata_dtacq',num2str(dtacq_num(2)),'_shot',num2str(shot(2)),'_tfshot',num2str(tfshot(2)),'.mat');
load(filename2,'rawdata');%a039
rawdata2=rawdata;
clear rawdata

%正しくデータ取得できていない場合はreturn
if numel(rawdata1)< 500||numel(rawdata2)<500
    return
end

%較正係数のバージョンを日付で判別
sheets1 = sheetnames('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\coeff125ch.xlsx');
sheets1 = str2double(sheets1);
sheet_date1=max(sheets1(sheets1<=date));
sheets2 = sheetnames('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\coeff200ch.xlsx');
sheets2 = str2double(sheets2);
sheet_date2=max(sheets2(sheets2<=date));

C1 = readmatrix('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\coeff125ch.xlsx','Sheet',num2str(sheet_date1));
C2 = readmatrix('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\coeff200ch.xlsx','Sheet',num2str(sheet_date2));

%a039
ok2 = logical(C2(:,14));
P2=C2(:,13);
coeff2=C2(:,12);
zpos2=C2(:,9);
rpos2=C2(:,10);
ch2=C2(:,7);

b2=rawdata2.*coeff2';%較正係数RC/NS
b2=b2.*P2';%極性揃え
b2=smoothdata(b2,1);

%デジタイザchからプローブ通し番号順への変換
bz2=zeros(1000,100);
bt2=bz2;
ok_bz2=false(100,1);
ok_bt2=ok_bz2;
zpos_bz2=zeros(100,1);
rpos_bz2=zpos_bz2;
zpos_bt2=zpos_bz2;
rpos_bt2=zpos_bz2;

for i=1:192
    if rem(ch2(i),2)==1
        bz2(:,ceil(ch2(i)/2))=b2(:,i);
        ok_bz2(ceil(ch2(i)/2))=ok2(i);
        zpos_bz2(ceil(ch2(i)/2))=zpos2(i);
        rpos_bz2(ceil(ch2(i)/2))=rpos2(i);
    elseif rem(ch2(i),2)==0
        bt2(:,ch2(i)/2)=b2(:,i);
        ok_bt2(ceil(ch2(i)/2))=ok2(i);
        zpos_bt2(ceil(ch2(i)/2))=zpos2(i);
        rpos_bt2(ceil(ch2(i)/2))=rpos2(i);
    end
end
[bz2, ok_bz2, ok_bz_plot2] = ng_replace(bz2, ok_bz2, sheet_date2);
% ok_bz2([48 58 49 59])=false;

%中心領域4+2本のみ
prange=21:80;%31:70;
bz2=bz2(:,prange);%time1000×ch225
zpos_bz2=zpos_bz2(prange);
rpos_bz2=rpos_bz2(prange);
ok_bz2=ok_bz2(prange);

%a038
ok1 = logical(C1(:,14));
P1=C1(:,13);
coeff1=C1(:,12);
zpos1=C1(:,9);
rpos1=C1(:,10);
ch1=C1(:,7);
coeff039= readmatrix('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\coeff125ch.xlsx','Sheet','a039coeff');
coeff039=[coeff039(:,1);coeff039(:,2);coeff039(:,3);coeff039(:,4);coeff039(:,5)];
coeff_v= readmatrix('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\coeff125ch.xlsx','Sheet','coeff_vacuum');
% p_ch= readmatrix('coeff125ch.xlsx','Sheet','p_ch');

b1=rawdata1.*coeff1'./1.1976;%較正係数RC/NS
b1=b1.*P1';%極性揃え
b1=double(b1);

%デジタイザchからプローブ通し番号順への変換
bz1=zeros(1000,126);
ok_bz1=true(100,1);
rpos_bz1=zeros(126,1);
zpos_bz1=rpos_bz1;

for i=1:128
    if ch1(i)>0
        bz1(:,ch1(i))=b1(:,i);
        ok_bz1(ch1(i))=ok1(i);
        rpos_bz1(ch1(i))=rpos1(i);
        zpos_bz1(ch1(i))=zpos1(i);
    end
end
bz1(:,63)=[];
ok_bz1(63)=[];
rpos_bz1(63)=[];
zpos_bz1(63)=[];

bz1=smoothdata(bz1,1);

ok_bz1([64 65 76])=false;
% ok_bz1([93 94])=false;
ok_bz1([19 31 38 83])=false;
% ok_bz1([40:43])=false;

bz1=bz1./coeff039';
%r<=0.1m以下は使わない
ok_bz1([1:3 26:28 51:53 76:78 101:103])=false;
rpos_bz1=rpos_bz1+0.008;

% bz1=bz1.*coeff_v';

% for i=1:125
%     bz1(:,i)=lowpass(bz1(:,i),0.5e5,1e6);
% end
% 

% ok_bz1([5 54 31 8 61 37 62 13 38 16 40 20 45 47 23 48 25 28 29 30])=false;
% ok_bz1([1:7 13:16 21:24 26:34 45:48 55:57 60 41:44])=false;
% bz1(1,:)=[];
% bz1=vertcat(bz1,zeros(1,125));


%データ統合
bz=[bz1 bz2];%time1000×ch(125+40)
zpos_bz=[zpos_bz1; zpos_bz2];
rpos_bz=[rpos_bz1; rpos_bz2];
ok_bz=[ok_bz1;ok_bz2];
% bz=smoothdata(bz,1);

% EF coil パラメータ 
r_EF   = 0.5 ;
n_EF   = 234. ;
if date<221119
    z1_EF   = 0.875;%0.68;
    z2_EF   = -0.830;%-0.68;
else
    z1_EF   = 0.78;
    z2_EF   = -0.78;
end

%125+40ch
[zq,rq]=meshgrid(linspace(min(zpos_bz),max(zpos_bz),n),linspace(min(rpos_bz),max(rpos_bz),n));
grid2D=struct('zq',zq,'rq',rq);
clear zq rq

[Bz_EF,~] = B_EF(z1_EF,z2_EF,r_EF,i_EF,n_EF,grid2D.rq,grid2D.zq,false);
% clear EF r_EF n_EF i_EF z_EF

data2D=struct('psi',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),'Bz',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),'Br',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),'Jt',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),'Et',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),'trange',trange);

for i=1:size(trange,2)
    t=trange(i);
    %%Bzの二次元補間(線形)
    vq = bz_rbfinterp(rpos_bz, zpos_bz, grid2D, bz, ok_bz, t);
    vq=imgaussfilt(vq);
    B_z = -Bz_EF+vq;
    %%PSI計算
    data2D.psi(:,:,i) = cumtrapz(grid2D.rq(:,1),2.*pi.*grid2D.rq(:,1).*B_z,1);
    %このままだと1/2πrが計算されてないので
    [data2D.Br(:,:,i),data2D.Bz(:,:,i)]=gradient(data2D.psi(:,:,i),grid2D.zq(1,:),grid2D.rq(:,1)) ;
    data2D.Br(:,:,i)=-data2D.Br(:,:,i)./(2.*pi.*grid2D.rq);
    data2D.Bz(:,:,i)=data2D.Bz(:,:,i)./(2.*pi.*grid2D.rq);
    data2D.Jt(:,:,i)= curl(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),data2D.Br(:,:,i))./(4*pi*1e-7);
end
data2D.Et=diff(data2D.psi,1,3).*1e+6; 
%diffは単なる差分なので時間方向のsizeが1小さくなる %ステップサイズは1us
data2D.Et=-1.*data2D.Et./(2.*pi.*grid2D.rq);

if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end
ok_z = zpos_bz(ok_bz); %z方向の生きているチャンネル
ok_r = rpos_bz(ok_bz); %r方向の生きているチャンネル

figure('Position', [0 0 1500 1500],'visible','on');
start=30;
%  t_start=470+start;
 for m=1:10 %図示する時間
     i=start+m; %end
     t=trange(i);
     subplot(2,5,m)
    contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),30,'LineStyle','none')
%     contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),10,'LineStyle','none')
    colormap(jet)
    axis image
    axis tight manual
%     caxis([-2.5*1e+6,2.5*1e+6]) %カラーバーの軸の範囲
    caxis([-0.04,0.04])
    colorbar('Location','eastoutside')
    hold on
%     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),30,'black')
    contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),[-20e-3:0.2e-3:40e-3],'black')
    plot(ok_z,ok_r,"k.",'MarkerSize', 7)%測定位置
    hold off
    title(string(t)+' us')
    xlabel('z [m]')
    ylabel('r [m]')
 end

 %a038:125ch
[zq1,rq1]=meshgrid(linspace(min(zpos_bz1),max(zpos_bz1),n),linspace(min(rpos_bz1),max(rpos_bz1),n));
grid2D1=struct('zq',zq1,'rq',rq1);
clear zq1 rq1

[Bz_EF1,~] = B_EF(z1_EF,z2_EF,r_EF,i_EF,n_EF,grid2D1.rq,grid2D1.zq,false);

data2D1=struct('psi',zeros(size(grid2D1.rq,1),size(grid2D1.rq,2),size(trange,2)),'Bz',zeros(size(grid2D1.rq,1),size(grid2D1.rq,2),size(trange,2)),'Br',zeros(size(grid2D1.rq,1),size(grid2D1.rq,2),size(trange,2)),'Jt',zeros(size(grid2D1.rq,1),size(grid2D1.rq,2),size(trange,2)),'Et',zeros(size(grid2D1.rq,1),size(grid2D1.rq,2),size(trange,2)),'trange',trange);

for i=1:size(trange,2)
    t=trange(i);
    %%Bzの二次元補間(線形)
    vq1 = bz_rbfinterp(rpos_bz1, zpos_bz1, grid2D1, bz1, ok_bz1, t);
    B_z1 = -Bz_EF1+vq1;
    %%PSI計算
    data2D1.psi(:,:,i) = cumtrapz(grid2D1.rq(:,1),2.*pi.*grid2D1.rq(:,1).*B_z1,1);
    %このままだと1/2πrが計算されてないので
    [data2D1.Br(:,:,i),data2D1.Bz(:,:,i)]=gradient(data2D1.psi(:,:,i),grid2D1.zq(1,:),grid2D1.rq(:,1)) ;
    data2D1.Br(:,:,i)=-data2D1.Br(:,:,i)./(2.*pi.*grid2D1.rq);
    data2D1.Bz(:,:,i)=data2D1.Bz(:,:,i)./(2.*pi.*grid2D1.rq);
    data2D1.Jt(:,:,i)= curl(grid2D1.zq(1,:),grid2D1.rq(:,1),data2D1.Bz(:,:,i),data2D1.Br(:,:,i))./(4*pi*1e-7);
end
data2D1.Et=diff(data2D1.psi,1,3).*1e+6; 
%diffは単なる差分なので時間方向のsizeが1小さくなる %ステップサイズは1us
data2D1.Et=-1.*data2D1.Et./(2.*pi.*grid2D1.rq);

if isstruct(grid2D1)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end

%a039:40ch
[zq2,rq2]=meshgrid(linspace(min(zpos_bz2),max(zpos_bz2),n),linspace(min(rpos_bz2),max(rpos_bz2),n));
grid2D2=struct('zq',zq2,'rq',rq2);
clear zq2 rq2

[Bz_EF2,~] = B_EF(z1_EF,z2_EF,r_EF,i_EF,n_EF,grid2D2.rq,grid2D2.zq,false);

data2D2=struct('psi',zeros(size(grid2D2.rq,1),size(grid2D2.rq,2),size(trange,2)),'Bz',zeros(size(grid2D2.rq,1),size(grid2D2.rq,2),size(trange,2)),'Br',zeros(size(grid2D2.rq,1),size(grid2D2.rq,2),size(trange,2)),'Jt',zeros(size(grid2D2.rq,1),size(grid2D2.rq,2),size(trange,2)),'Et',zeros(size(grid2D2.rq,1),size(grid2D2.rq,2),size(trange,2)),'trange',trange);

for i=1:size(trange,2)
    t=trange(i);
    %%Bzの二次元補間(線形)
    vq2 = bz_rbfinterp(rpos_bz2, zpos_bz2, grid2D2, bz2, ok_bz2, t);
    B_z2 = -Bz_EF2+vq2;
    %%PSI計算
    data2D2.psi(:,:,i) = cumtrapz(grid2D2.rq(:,1),2.*pi.*grid2D2.rq(:,1).*B_z2,1);
    %このままだと1/2πrが計算されてないので
    [data2D2.Br(:,:,i),data2D2.Bz(:,:,i)]=gradient(data2D2.psi(:,:,i),grid2D2.zq(1,:),grid2D2.rq(:,1)) ;
    data2D2.Br(:,:,i)=-data2D2.Br(:,:,i)./(2.*pi.*grid2D2.rq);
    data2D2.Bz(:,:,i)=data2D2.Bz(:,:,i)./(2.*pi.*grid2D2.rq);
    data2D2.Jt(:,:,i)= curl(grid2D2.zq(1,:),grid2D2.rq(:,1),data2D2.Bz(:,:,i),data2D2.Br(:,:,i))./(4*pi*1e-7);
end
data2D2.Et=diff(data2D2.psi,1,3).*1e+6; 
%diffは単なる差分なので時間方向のsizeが1小さくなる %ステップサイズは1us
data2D2.Et=-1.*data2D2.Et./(2.*pi.*grid2D2.rq);

if isstruct(grid2D2)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end

figure('Position', [0 0 1500 1500],'visible','on');
start=30;
 for m=1:5 %図示する時間
     i=start+m; %end
     t=trange(i);
     subplot(2,5,m.*2-1)
    contourf(grid2D1.zq(1,:),grid2D1.rq(:,1),data2D1.Bz(:,:,i),30,'LineStyle','none')
%     contourf(grid2D1.zq(1,:),grid2D1.rq(:,1),data2D1.psi(:,:,i),30,'LineStyle','none')
%     contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),10,'LineStyle','none')
    colormap(jet)
    axis image
    axis tight manual
%     caxis([-2.5*1e+6,2.5*1e+6]) %カラーバーの軸の範囲
%     caxis([0,5e-3])
        caxis([-0.04,0.04])
%     colorbar('Location','eastoutside')
    hold on
%     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),30,'black')
    contour(grid2D1.zq(1,:),grid2D1.rq(:,1),squeeze(data2D1.psi(:,:,i)),[-20e-3:0.2e-3:40e-3],'black')
%     plot(ok_z,ok_r,"k.",'MarkerSize', 7)%測定位置
    hold off
    xlim([-0.085 0.085])%xlim([-0.0525 0.0525])
    ylim([0.06 0.33])
    title(string(t)+' us')
    
    subplot(2,5,m.*2)
    contourf(grid2D2.zq(1,:),grid2D2.rq(:,1),data2D2.Bz(:,:,i),30,'LineStyle','none')
%     contourf(grid2D2.zq(1,:),grid2D2.rq(:,1),data2D2.psi(:,:,i),30,'LineStyle','none')
%     contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),10,'LineStyle','none')
    colormap(jet)
    axis image
    axis tight manual
%     caxis([-2.5*1e+6,2.5*1e+6]) %jt
%     caxis([0,5e-3])%psi
    caxis([-0.04,0.04])%Bz
%     colorbar('Location','eastoutside')
    hold on
%     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),30,'black')
    contour(grid2D2.zq(1,:),grid2D2.rq(:,1),squeeze(data2D2.psi(:,:,i)),[-20e-3:0.2e-3:40e-3],'black')
%     plot(ok_z,ok_r,"k.",'MarkerSize', 7)%測定位置
    hold off
    xlim([-0.085 0.085])%xlim([-0.0525 0.0525])
    ylim([0.06 0.33])
    title(string(t)+' us')
    
 end

end
