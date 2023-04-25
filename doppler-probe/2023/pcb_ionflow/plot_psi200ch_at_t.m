%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%プロット時間を指定して磁気面を1枚プロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = plot_psi200ch_at_t(time,date,dtacq,pathname,mesh_rz,expval,trange,cmap)
%プロット時刻/実験日/dtacq番号(39)/dtacqshot番号/dtacqtfshot番号/pathname/rz方向メッシュ数/EF電流/計算時間範囲/カラーマップ(TorF)

filename=strcat(pathname.rawdata,'/',num2str(date),'/rawdata_dtacq',num2str(dtacq.num),'_shot',num2str(dtacq.shot),'_tfshot',num2str(dtacq.tfshot),'.mat');
if exist(filename,"file")==0
    warning(strcat(filename,' does not exist.'));
    return
end
load(filename,'rawdata');%1000×192

%正しくデータ取得できていない場合はreturn
if numel(rawdata)< 500
    disp('data incomplete/corrupted')
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
ok_bt([4 5 6 7 8 9 10 15 21 27 30 42 43 49 53 69 84 87 92 94 95 96 97 98 99 100]) = false;

[zq,rq]=meshgrid(linspace(min(zpos_bz),max(zpos_bz),mesh_rz),linspace(min(rpos_bz),max(rpos_bz),mesh_rz));

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
[Bz_EF,~] = B_EF(z1_EF,z2_EF,r_EF,expval.EF,n_EF,grid2D.rq,grid2D.zq,false);
clear EF r_EF n_EF expval.EF z_EF

data2D=struct('psi',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Bz',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Bt',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Br',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Jt',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
    'Et',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),'trange',trange);

for i=1:size(trange,2)
    t=trange(i);
    %%Bzの二次元補間(線形fit)
    vq = bz_rbfinterp(rpos_bz, zpos_bz, grid2D, bz, ok_bz, t);
    B_z = -Bz_EF+vq;
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

figure('Position', [0 0 1500 1500],'visible','on');
%  t_start=470+start;
i=(time-trange(1)+1);
t=trange(i);
if cmap
    %     contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),30,'LineStyle','none')
    %     contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.psi(:,:,i),40,'LineStyle','none')
    contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt(:,:,i),-100e-3:0.5e-3:100e-3,'LineStyle','none')
    % contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),30,'LineStyle','none')
    %     contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Et(:,:,i),20,'LineStyle','none')
    colormap(jet)
    axis image
    axis tight manual
    %     caxis([-0.8*1e+6,0.8*1e+6]) %jt%カラーバーの軸の範囲
    %     caxis([-0.01,0.01])%Bz
    caxis([-0.1,0.1])%Bt
    %     caxis([-5e-3,5e-3])%psi
    %     caxis([-500,400])%Et
    %     colorbar('Location','eastoutside')
    %カラーバーのラベル付け
    %     c = colorbar;
    %     c.Label.String = 'Jt [A/m^{2}]';
    hold on
end
%     plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
%     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')
%     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')
contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),-20e-3:0.2e-3:40e-3,'black','LineWidth',1)
%     plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"bo")
%     plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"bx")
hold on
plot(ok_z,ok_r,"k.",'MarkerSize', 6)%測定位置
hold on
title(string(t)+' us')
xlabel('z [m]')
ylabel('r [m]')
ax = gca;
ax.FontSize = 14;

% figure('Position', [0 0 1500 1500],'visible','on');
% for m=1:16 %図示する時間
%     i=(start-trange(1)+1)+(m-1)*dt; %end
%     t=trange(i);
%     subplot(4,4,m)
%     plot(grid2D.zq(1,:),data2D.Bt(25,:,i))
%     title(string(t)+' us')
%     ylim([-0.03 0.03])
% end
% saveas(gcf,strcat('/Users/yunhancai/Downloads/files/a039_',num2str(shot)))
% saveas(gcf,strcat('/Users/yunhancai/Downloads/a039_',num2str(shot),'.png'))
% close
end
