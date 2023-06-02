%%%%%%%%%%%%%%%%%%%%%%%%
%200ch用新規pcbプローブの電流シート計算
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%ここが各PCのパス
%【※コードを使用する前に】環境変数を設定しておくか、matlab内のコマンドからsetenv('パス名','アドレス')で指定してから動かす
pathname.rawdata=getenv('rawdata_path');%dtacqのrawdataの保管場所

%%%%実験オペレーションの取得
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
pat=230127;
T=searchlog(T,node,pat);
%IDXlist=[4:6 8:11 13 15:19 21:23 24:30 33:37 39:40 42:51 53:59 61:63 65:69 71:74];
IDXlist = 63; 
date=pat;
n_data=numel(IDXlist);%計測データ数
shotlist=T.a039(IDXlist);
tfshotlist=T.a039_TF(IDXlist);
EFlist=T.EF_A_(IDXlist);
TFlist=T.TF_kV_(IDXlist);
dtacqlist=39.*ones(n_data,1);

% % %直接入力の場合
% dtacqlist=39;
% shotlist=413;%240;%【input】dtacqの保存番号
% tfshotlist=411;%0;
% date = 230119;%【input】計測日
% n_data=numel(shotlist);%計測データ数
% EFlist = 150;%150;%【input】EF電流
% TFlist=4;

trange=450:510;%【input】計算時間範囲
n=50; %【input】rz方向のメッシュ数

for i=1:n_data
    dtacq_num=dtacqlist(i);
    shot=shotlist(i);
    tfshot=tfshotlist(i);
    i_EF=EFlist(i);
    TF=TFlist(i);
    [p,max_slope,min_slope]=calc_tilt(date, dtacq_num, shot, tfshot, pathname,n,i_EF,trange,TF); 
end

%%%%%%%%%%%%%%%%%%%%%%%%
%以下、local関数
%%%%%%%%%%%%%%%%%%%%%%%%

function [p,max_slope,min_slope]=calc_tilt(date, dtacq_num, shot, tfshot, pathname, n,i_EF,trange,TF)
filename=strcat(pathname.rawdata,'rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat');
if exist(filename,"file")==0
    return
end
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
[bz, ok_bz, ok_bz_plot] = ng_replace(bz, ok_bz, sheet_date);
% ok_bz_plot=ok_bz;

%中心領域4+2本のみ
prange=21:80;%31:70;
bz=bz(:,prange);%time1000×ch225
zpos_bz=zpos_bz(prange);
rpos_bz=rpos_bz(prange);
ok_bz=ok_bz(prange);
ok_bz_plot=ok_bz_plot(prange);

[zq,rq]=meshgrid(linspace(min(zpos_bz),max(zpos_bz),n),linspace(min(rpos_bz),max(rpos_bz),n));
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
clear EF r_EF n_EF z_EF

data2D=struct('psi',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),'Bz',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),'Br',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),'Jt',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),'Et',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),'trange',trange);

for i=1:size(trange,2)
    t=trange(i);
    %%Bzの二次元補間(線形fit)
    vq =bz_rbfinterp(rpos_bz, zpos_bz, grid2D, bz, ok_bz, t);
    B_z = -Bz_EF+vq;
    %%PSI計算
    data2D.psi(:,:,i) = cumtrapz(grid2D.rq(:,1),2*pi*B_z.*grid2D.rq(:,1),1);
    %このままだと1/2πrが計算されてないので
    [data2D.Br(:,:,i),data2D.Bz(:,:,i)]=gradient(data2D.psi(:,:,i),grid2D.zq(1,:),grid2D.rq(:,1)) ;
    data2D.Br(:,:,i)=-data2D.Br(:,:,i)./(2.*pi.*grid2D.rq);
    data2D.Bz(:,:,i)=data2D.Bz(:,:,i)./(2.*pi.*grid2D.rq);
    data2D.Jt(:,:,i)= curl(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),data2D.Br(:,:,i))./(4*pi*1e-7);
end
data2D.Et=diff(data2D.psi,1,3).*1e+6; 
%diffは単なる差分なので時間方向のsizeが1小さくなる %ステップサイズは1us
data2D.Et=-1.*data2D.Et./(2.*pi.*grid2D.rq);

% ok_z = zpos_bz(ok_bz_plot); %z方向の生きているチャンネル
% ok_r = rpos_bz(ok_bz_plot); %r方向の生きているチャンネル

if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end

[psimid,mid]=min(data2D.psi,[],2);%磁気中性面
[opoint,~]=islocalmin(psimid,1);
[xpoint,~]=islocalmax(psimid,1);
r_index=1:n;

frame=numel(trange);

row = 2;
column =ceil(frame/row);
%X点での勾配
p=zeros(2,frame);
% figure
% f1=figure;
% f1.WindowState='maximized';
for i=1:frame
%     r_xind=ceil(n/2);
%     z_xind=ceil(n/2);
    if sum(xpoint(:,:,i))>0   
        r_xind=r_index(xpoint(:,:,i));
        z_xind=mid(xpoint(:,:,i),:,i);
        if numel(r_xind)>1
            [~,I]=min(psimid(r_xind));
            r_xind=r_xind(I);
            z_xind=z_xind(I);
        end
        
        r_s=max(2,r_xind-5);
        r_n=min(n,r_xind+5);
        p(:,i)=polyfit(grid2D.rq(r_s:r_n,1),grid2D.zq(1,squeeze(mid(r_s:r_n,:,i))),1);
        tilt=polyval(p(:,i),grid2D.rq(r_s:r_n,1));

%         subplot(row,column,i)
%         contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),30,'LineStyle','none');
%         colormap(jet)
%         caxis([-0.7*1e+6,0.7*1e+6]) 
%         axis image
%         axis tight manual
%         hold on
%         contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),[-20e-3:0.2e-3:40e-3],'black','LineWidth',0.3)
% %         colorbar('Location','eastoutside')
%         plot(tilt,grid2D.rq(r_s:r_n,1),'r-','LineWidth',2)
%         plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"bx",'MarkerSize',8)
%         plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1),'b-','LineWidth',1)
%         hold off
%         alpha(0.3)
%         ylim([0.1 grid2D.rq(end,1)])
%         title(string(trange(1)+i-1)+' us')

        if z_xind==1 || z_xind==n
            p(:,i)=NaN;
        end
    else
        p(:,i)=NaN;
    end   
end

figure
plot(trange,p(1,:).*180./pi,'o-','LineWidth',1)
hold on
yline(0,'k--')
xlim([450 500])
ylim([-15 15])
% ylim([-0.55 0.1].*180./pi)
xlabel('Time [us]')
ylabel('slope near X point')
% saveas(gcf,strcat(pathname.save,'\IDX',num2str(IDX),'slope.png'))

max_slope=max(abs(p(1,:)));
min_slope=min(p(1,:));

saveas(gcf,strcat('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\230123\sheet\',num2str(shot),'_sheetangle_TF',num2str(TF),'kV.png'))
save(strcat('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\230123\sheet\sheetangle_',num2str(shot),'.mat'),'date','i_EF','n','shot','TF','tfshot','trange','p','min_slope','max_slope')
close

% %各時間、各列(z)ごとのpsiの最大値
% [max_psi,max_psi_r]=max(data2D.psi,[],1);
% max_psi=squeeze(max_psi);
% psi_pr=zeros(3,frame);
% xJt=zeros(1,frame);
% xEt=xJt;
% xpos=zeros(2,frame);%X点の座標(z,r)
% 
% for i=1:frame
%     r_ind=max_psi_r(:,:,i);
%     max_psi_ind=find(islocalmax(smooth(max_psi(:,i)),'MaxNumExtrema', 2));
%     if numel(max_psi(min(max_psi_ind),i))==0
%     psi_pr(1,i)=NaN;%psi_pr(1,i)=max_psi(1,i);
%     else
%     psi_pr(1,i)=max_psi(min(max_psi_ind),i);
%     end
%     if numel(max_psi(max(max_psi_ind),i))==0
%     psi_pr(2,i)=NaN;%psi_pr(2,i)=max_psi(n,i);
%     else
%     psi_pr(2,i)=max_psi(max(max_psi_ind),i);
%     end
%     if numel(find(islocalmin(smooth(max_psi(:,i)),'MaxNumExtrema', 1)))==0
%         psi_pr(3,i)=NaN;
%         xJt(1,i)=NaN;
%         xEt(1,i)=NaN;
%         xpos(:,i)=NaN;
%     else
%         min_psi_ind=islocalmin(smooth(max_psi(:,i)),'MaxNumExtrema', 1);
%         xr=r_ind(min_psi_ind);
%         if xr==1 || xr==n %r両端の場合は検知しない
%             psi_pr(3,i)=NaN;
%             xJt(1,i)=NaN;
%             xEt(1,i)=NaN;
%             xpos(:,i)=NaN;
%         else
%         psi_pr(3,i)=max_psi(min_psi_ind,i);
%         xJt(1,i)=data2D.Jt(xr,min_psi_ind,i);
%         xEt(1,i)=data2D.Et(xr,min_psi_ind,i);
%         xpos(1,i)=grid2D.zq(1,min_psi_ind);
%         xpos(2,i)=grid2D.rq(xr,1);
%         end
%     end
%     if xJt(1,i)>-1e3
%         xJt(1,i)=NaN;
%     end
% 
% %     if max_psi(1,i)==max(max_psi(1:end/2,i))
% %         psi_pr(1,i)=NaN; %max_psi(1,i);
% %     end    
% %     if max_psi(end,i)==max(max_psi(end/2:end,i))
% %         psi_pr(2,i)=NaN; %max_psi(end,i);
% %     end
% end
% fitrate=psi_pr(3,:)./min(psi_pr(1:2,:),[],1);
% xeta=xEt(1,:)./xJt(1,:);
% xeta(xeta>0.01)=NaN;
% 
% %合体率とX点Jt,Et,etaのplot
% figure('Position',[0,0,300,500])
% 
% subplot(2,1,1)
% yyaxis left
% xJt_plot = rmmissing(xJt);
% trange_plot=trange(not(isnan(xJt)));
% plot(trange_plot,-1.*xJt_plot,'b-+','LineWidth',1)
% ylabel('Jt [A/m^{2}]')
% xlabel('Time [us]')
% xlim([450 500])
% ylim([0 4e5])
% yyaxis right
% plot(trange,-1.*xEt,'r-+','LineWidth',1)
% ylabel('Et [V/m]')
% xlim([450 500])
% ylim([0 300])
% ha2 = gca;
% ha2.LineWidth = 1;
% ha2.FontSize=10;
% 
% subplot(2,1,2)
% xeta_plot = rmmissing(xeta);
% trange_plot=trange(not(isnan(xeta)));
% plot(trange_plot,xeta_plot,'k-+','LineWidth',1)
% ylabel('η [Ω m]')
% xlabel('Time [us]')
% xlim([450 500])
% ylim([0 inf])
% ha1 = gca;
% ha1.LineWidth = 1;
% ha1.FontSize=10;
% 
% % saveas(gcf,strcat('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\230123\xpoint\',num2str(shot),'_xpoint_TF',num2str(TF),'kV.png'))
% % save(strcat('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\230123\xpoint\xpoint_',num2str(shot),'.mat'),'date','i_EF','fitrate','n','psi_pr','shot','TF','tfshot','trange','xEt','xeta','xJt','xpos')
% % close
% 
% %X点と磁気面の重ね合わせ
% figure('Position', [0 0 1500 1500],'visible','on');
% start=15;
% %  t_start=470+start;
%  for m=1:10 %図示する時間
%      i=start+m.*2; %end
%      t=trange(i);
%      subplot(2,5,m)
% %     contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),30,'LineStyle','none')
% %     contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.psi(:,:,i),40,'LineStyle','none')
%     contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),30,'LineStyle','none')
% %     contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Et(:,:,i),20,'LineStyle','none')
%     colormap(jet)
%     axis image
%     axis tight manual
%     caxis([-0.7*1e+6,0.7*1e+6]) %jt%カラーバーの軸の範囲
% %     caxis([-0.07,0.07])%Bz
% %     caxis([-5e-3,5e-3])%psi
% %     caxis([-500,400])%Et
% %     colorbar('Location','eastoutside')
%     %カラーバーのラベル付け
% %     c = colorbar;
% %     c.Label.String = 'Jt [A/m^{2}]';
%     hold on
%     plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
% % contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')
% % contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),40,'black')
%     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),[-20e-3:0.2e-3:40e-3],'black')
%     plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"ko")
%     plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"kx")
% %     plot(ok_z,ok_r,"k.",'MarkerSize', 6)%測定位置
%     plot(xpos(1,i),xpos(2,i),'bx')%X点
%     hold off
%     title(string(t)+' us')
% %     xlabel('z [m]')
% %     ylabel('r [m]')
%  end
% % %  saveas(gcf,strcat('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\230123\psi-jt\',num2str(shot),'_jt_TF',num2str(TF),'kV_t466us.png'))
% % %  close
end
