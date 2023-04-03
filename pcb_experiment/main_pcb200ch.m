%%%%%%%%%%%%%%%%%%%%%%%%
%200ch用新規pcbプローブのみでの磁気面（Bz）
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
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
pat=230127;
T=searchlog(T,node,pat);
IDXlist=63;%[7:9 12];%230111;[6 10 13:16 18:20 22 23 25:32 34:44 47:49 58 60 63 68:81];%230128%[4:6 8:11 13 15:19 21:23 24:30 33:37 39:40 42:51 53:59 61:63 65:69 71:74];%230119
date=pat;
n_data=numel(IDXlist);%計測データ数
shotlist=T.a039(IDXlist);
tfshotlist=T.a039_TF(IDXlist);
EFlist=T.EF_A_(IDXlist);
TFlist=T.TF_kV_(IDXlist);
dtacqlist=39.*ones(n_data,1);

% % %直接入力の場合
% dtacqlist=39;
% shotlist=650;%240;%【input】dtacqの保存番号
% tfshotlist=584;%0;
% date = 230128;%【input】計測日
% n_data=numel(shotlist);%計測データ数
% EFlist = 150;%150;%【input】EF電流
% TFlist=4;

trange=430:500;%【input】計算時間範囲
n=100; %【input】rz方向のメッシュ数

for i=1:n_data
    dtacq_num=dtacqlist(i);
    shot=shotlist(i);
    tfshot=tfshotlist(i);
    i_EF=EFlist(i);
%     TF=TFlist(i);
    plot_psi200ch(date, dtacq_num, shot, tfshot, pathname,n,i_EF,trange); 
end


%%%%%%%%%%%%%%%%%%%%%%%%
%以下、local関数
%%%%%%%%%%%%%%%%%%%%%%%%

function plot_psi200ch(date, dtacq_num, shot, tfshot, pathname, n,i_EF,trange)
filename=strcat(pathname.rawdata,'\rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat');
load(filename,'rawdata');%1000×192

%正しくデータ取得できていない場合はreturn
if numel(rawdata)< 500
    return
end

%較正係数のバージョンを日付で判別
sheets = sheetnames('C:\Users\uswk0\OneDrive\デスクトップ\Github\test-open\pcb_experiment\coeff200ch.xlsx');
sheets = str2double(sheets);
sheet_date=max(sheets(sheets<=date));

C = readmatrix('C:\Users\uswk0\OneDrive\デスクトップ\Github\test-open\pcb_experiment\coeff200ch.xlsx','Sheet',num2str(sheet_date));
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
ok_bz=false(100,1);
ok_bt=ok_bz;
zpos_bz=zeros(100,1);
rpos_bz=zpos_bz;
zpos_bt=zpos_bz;
rpos_bt=zpos_bz;

for i=1:192
    if rem(ch(i),2)==1
        bz(:,ceil(ch(i)/2))=b(:,i);
        ok_bz(ceil(ch(i)/2))=ok(i);
        zpos_bz(ceil(ch(i)/2))=zpos(i);
        rpos_bz(ceil(ch(i)/2))=rpos(i);
    elseif rem(ch(i),2)==0
        bt(:,ch(i)/2)=b(:,i);
        ok_bt(ceil(ch(i)/2))=ok(i);
        zpos_bt(ceil(ch(i)/2))=zpos(i);
        rpos_bt(ceil(ch(i)/2))=rpos(i);
    end
end
% ok_bz([9 10 96 17])=false;
% ok_bz([52 54 43])=false;
% ok_bz(57)=true;
[bz, ok_bz, ok_bz_plot] = ng_replace(bz, ok_bz, sheet_date);
% ok_bz_plot=ok_bz;
% ok_bz([48 58 49 59])=false;

% ok_bz([44 45])=false;

%中心領域4+2本のみ
prange=21:80;%31:70;
bz=bz(:,prange);%time1000×ch225
zpos_bz=zpos_bz(prange);
rpos_bz=rpos_bz(prange);
ok_bz=ok_bz(prange);
ok_bz_plot=ok_bz_plot(prange);


[zq,rq]=meshgrid(linspace(min(zpos_bz),max(zpos_bz),n),linspace(min(rpos_bz),max(rpos_bz),n));
% [zq,rq]=meshgrid(linspace(-0.0525,0.0525,51),linspace(0,0.25,51));
grid2D=struct('zq',zq,'rq',rq);
clear zq rq

%data2Dcalc.m
r_EF   = 0.5 ;
n_EF   = 234. ;
%i_EF    =i_EF;

if date<221119
    z1_EF   = 0.875;%0.68;
    z2_EF   = -0.830;%-0.68;
else
    z1_EF   = 0.78;
    z2_EF   = -0.78;
end
[Bz_EF,~] = B_EF(z1_EF,z2_EF,r_EF,i_EF,n_EF,grid2D.rq,grid2D.zq,false);
clear EF r_EF n_EF i_EF z_EF

data2D=struct('psi',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)), ...
    'Bz',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)), ...
    'Bt',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Br',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)), ...
    'Jt',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)), ...
    'Et',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)), ...
    'trange',trange);

for i=1:size(trange,2)
    t=trange(i);
    %%Bzの二次元補間(線形fit)
    vq =bz_rbfinterp(rpos_bz, zpos_bz, grid2D, bz, ok_bz, t);
    B_z = -Bz_EF+vq;
    %%Btの二次元補間(線形fit)
    B_t = bz_rbfinterp(rpos_bt, zpos_bt, grid2D, bt, ok_bt, t);
    %%PSI計算
    data2D.psi(:,:,i) = cumtrapz(grid2D.rq(:,1),2*pi*B_z.*grid2D.rq(:,1),1);
    %このままだと1/2πrが計算されてないので
    [data2D.Br(:,:,i),data2D.Bz(:,:,i)]=gradient(data2D.psi(:,:,i),grid2D.zq(1,:),grid2D.rq(:,1)) ;
    data2D.Br(:,:,i)=-data2D.Br(:,:,i)./(2.*pi.*grid2D.rq);
    data2D.Bz(:,:,i)=data2D.Bz(:,:,i)./(2.*pi.*grid2D.rq);
    data2D.Bt(:,:,i)=B_t;
    data2D.Jt(:,:,i)= curl(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),data2D.Br(:,:,i))./(4*pi*1e-7);
end
data2D.Et=diff(data2D.psi,1,3).*1e+6; 
%diffは単なる差分なので時間方向のsizeが1小さくなる %ステップサイズは1us
data2D.Et=data2D.Et./(2.*pi.*grid2D.rq);

ok_z = zpos_bz(ok_bz_plot); %z方向の生きているチャンネル
ok_r = rpos_bz(ok_bz_plot); %r方向の生きているチャンネル

if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end

% figure
% for i=1:50
%     minEt(i)=min(data2D.Et(:,:,i),[],'all');
% end
% plot(440:489,minEt)

%figure('Position', [0 0 1500 1500],'visible','on');
h=figure;
h.WindowState = 'maximized';

start=30;
%  t_start=470+start;
 for m=1:16 %図示する時間
     i=start+m.*2; %end
     t=trange(i);
     subplot(4,4,m)
%     contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),30,'LineStyle','none')
%     contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.psi(:,:,i),50,'LineStyle','none')
    contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),20,'LineStyle','none')
%     contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt(:,:,i),-80e-3:0.4e-3:80e-3,'LineStyle','none')
%     contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Et(:,:,i),20,'LineStyle','none')
    colormap(jet)
    axis image
    axis tight manual
   caxis([-1.6*1e+6,1.3*1e+6])%caxis([-1.4*1e+6,1.4*1e+6])%%caxis([-1.9*1e+6,1.3*1e+6])%%jt%カラーバーの軸の範囲
%     caxis([-0.05,0.05])%Bz
%     caxis([-8e-3,8e-3])%psi
%     caxis([-500,400])%Et
%     colorbar('Location','eastoutside')
    %カラーバーのラベル付け
%     c = colorbar;
%     c.Label.String = 'Jt [A/m^{2}]';
    hold on

% plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
% contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),50,'black')
% contour(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),15,'w-','LineWidth',0.01)
% contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),65,'black')
% contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),[-20e-3:0.1e-3:40e-3],'black')
    contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),[-20e-3:0.2e-3:40e-3],'black')
%     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),[0.36e-3 1.09e-3],'black')
% contour(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),[0.7e6 0.64e6 0.53e6],'w-','LineWidth',0.01)



%     plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"bo")
%     plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"bx")
    plot(ok_z,ok_r,"k.",'MarkerSize', 6)%測定位置
    hold off
    title(string(t)+' us')
     xlabel('z [m]')
     ylabel('r [m]')
 end
%  saveas(gcf,strcat(''))
%  close
end
