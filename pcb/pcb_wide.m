%電流シートの幅を求める

clearvars

pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス
pathname.fourier=getenv('fourier_path');%md0までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath
pathname.save='C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\過去データ解析\221229';%220926
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';
T=getTS6log(DOCID);% ログのテーブルを取得

%(a)No Cr
IDXlist=2911;%2911:2913;
tstart=472;
tend=484;

% %(b)Cr 454us
% IDXlist=2925:2927;
% tstart=477;
% tend=487;
 
% %(c)Cr 450us
% IDXlist=[2931 2933];
% tstart=482;
% tend=500;

tsize = tend-tstart+1;
wide = zeros(tsize,3); %電流シート幅

for i=1:numel(IDXlist)
    IDX=IDXlist(i);
    wide(:,i)=cal_diff(T, IDX, pathname, tstart, tend, tsize);
end

function [date, shot, TF_shot, offset_TF, i_EF, start, v_TF, d_tacq, d_tacqTF, n] = getinput(T,IDX)
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

n=50; %rz方向のメッシュ数
end

function wide=cal_diff(T, IDX, pathname, tstart, tend, tsize)
trange = 470:510; %Jtは差分で求めるため、欲しい値の一つ前からとる
row = 2;
column =ceil(tsize/row);
[date, ~, ~, ~, i_EF, ~, ~, d_tacq, d_tacqTF, n] = getinput(T,IDX);%Tのテーブルから入力のリストを出力
[grid2D, data2D] = pcbdata(date, d_tacq, d_tacqTF, trange, [], n,i_EF);

%tstart:tendの範囲以外は削除（補正や微分の関係上、見たい時間ぴったりの解析では不適切と思ったため）
data2D.psi(:,:,tend-trange(1)+2:numel(trange))=[];
data2D.psi(:,:,1:tstart-trange(1))=[];
data2D.Bz(:,:,tend-trange(1)+2:numel(trange))=[];
data2D.Bz(:,:,1:tstart-trange(1))=[];
data2D.Br(:,:,tend-trange(1)+2:numel(trange))=[];
data2D.Br(:,:,1:tstart-trange(1))=[];
data2D.Et(:,:,tend-trange(1)+2:numel(trange)-1)=[];%Etはdiffで求めているので時間一つ分少ない
data2D.Et(:,:,1:tstart-trange(1))=[];
data2D.Jt(:,:,tend-trange(1)+2:numel(trange))=[];
data2D.Jt(:,:,1:tstart-trange(1))=[];

Jtmin=min(data2D.Jt,[],[1 2]);

%%%midplaneとかO点、X点を探す
[psimid,mid]=min(data2D.psi,[],2);
[opoint,~]=islocalmin(psimid,1);
[xpoint,~]=islocalmax(psimid,1);
% [xp_psi,maxxp]=max(squeeze(psimid),[],1);
r_index=1:n;

dr=grid2D.rq(2,1)-grid2D.rq(1,1);
dz=grid2D.zq(1,2)-grid2D.zq(1,1);
r_value=grid2D.rq(:,1);

Q_n=zeros(n,1);%シートの接線と各z-lineの上側交点
Q_s=Q_n;%シートの接線と各z-lineの下側交点
Jt=zeros(n,tsize);%垂線上のJt
Bp=Jt;%垂線上でmidplaneに平行な磁場成分(Br)
Bv=Jt;%垂線上でmidplaneに垂直な磁場成分(Bz)
p=zeros(2,tsize);%midplaneの近似直線の係数（引数r）
v=zeros(2,tsize);%midplaneの近似直線の接線の係数（引数z）
wide=zeros(tsize,1);%シート幅
r_xind=zeros(1,tsize);%X-point(1点)のr座標index
z_xind=r_xind;%X-point(1点)のz座標index
ri=zeros(2,tsize);%シート幅端のr座標(左右で2個)
zi=zeros(2,tsize);%シート幅端のz座標(左右で2個)
sheetJt=zeros(tsize);

f1=figure;
f1.WindowState='maximized';
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
        r_xind(1,i)=r_xind0;
        z_xind(1,i)=z_xind0;
        clear z_xind0 r_xind0

        %(2)X点の上下5マスのmidplaneの線形近似→シートに平行・垂直な単位ベクトル
        r_s=max(2,r_xind(:,i)-5);
        r_n=min(n,r_xind(:,i)+5);
        p(:,i)=polyfit(grid2D.rq(r_s:r_n,1),grid2D.zq(1,squeeze(mid(r_s:r_n,:,i))),1);
        diff=polyval(p(:,i),grid2D.rq(r_s:r_n,1));%midplaneの近似直線のz座標(1×11)

        v(:,i)=[-p(1,i) grid2D.rq(r_xind(:,i),1)+p(1,i).*grid2D.zq(1,z_xind(:,i))];
        vertical=polyval(v(:,i),grid2D.zq(1,1:n));%垂線のr座標(1×n)

        ep=[p(1,i) 1]./(p(1,i)^2+1); %シートに平行な単位ベクトル
        ev=[1 v(1,i)]./(v(1,i)^2+1); %シートに垂直な単位ベクトル

        %(3)垂線と各z-lineの交点に最も近いメッシュ点の座標(grid2D.zq(1,j),Q_n(j)), (grid2D.zq(1,j),Q_s(j))
        for j=1:n
            if r_value(n)>vertical(j)%傾きが大きく、シート幅右端がrmaxの端から出た場合を省く
               Q_n(j)=min(r_value(r_value>vertical(j)));
               Q_s(j)=max(r_value(r_value<=vertical(j)));
            else
                Q_n(j)=NaN;
                Q_s(j)=NaN;
            end 
        end

        Q_nind=round((Q_n-grid2D.rq(1,1))./dr)+1;%Q_nのindex
        Q_sind=round((Q_s-grid2D.rq(1,1))./dr)+1;%Q_sのindex

        Q_nind(isnan(Q_nind)) = [];
        Q_sind(isnan(Q_sind)) = [];
        nq=numel(Q_nind);%表示されるverticalの要素

        %(4)Q_n, Q_sでの値から垂線と各z-lineの交点での値を近似
        Bp_n=diag(data2D.Bz(Q_nind,:,i)).*ep(1)+diag(data2D.Br(Q_nind,:,i)).*ep(2);
        Bv_n=diag(data2D.Bz(Q_nind,:,i)).*ev(1)+diag(data2D.Br(Q_nind,:,i)).*ev(2);
        Jt_n=diag(data2D.Jt(Q_nind,:,i));
        Bp_s=diag(data2D.Bz(Q_sind,:,i)).*ep(1)+diag(data2D.Br(Q_sind,:,i)).*ep(2);
        Bv_s=diag(data2D.Bz(Q_sind,:,i)).*ev(1)+diag(data2D.Br(Q_sind,:,i)).*ev(2);
        Jt_s=diag(data2D.Jt(Q_sind,:,i));

        %内分から求める場合
        di=[vertical'-Q_s Q_n-vertical'];%internally divide
        di=di(1:nq,:);
        Jt(1:nq,i)=dot(di, [Jt_n Jt_s],2)./dr;
        Bp(1:nq,i)=dot(di, [Bp_n Bp_s],2)./dr;
        Bv(1:nq,i)=dot(di, [Bv_n Bv_s],2)./dr;
        if nq<n
           Jt(nq+1:end,i)=NaN;
           Bp(nq+1:end,i)=NaN;
           Bv(nq+1:end,i)=NaN;
        end

        [~,min_ind]=min(Jt(:,i));%垂線上でのJt最小の点のz-index

%         %%(5-1:poly)垂線上のJtとmin(jt)*0.5の交点による半値全幅FWHM
%         % polyxpoly：Jtの各点を直線近似した場合の交点
%         [zi0,yi0,ii0]=polyxpoly(grid2D.zq(1,1:nq),Jt(1:nq,i),grid2D.zq(1,1:nq),ones(nq,1).*min(Jt(:,i)).*0.5);
%         if numel(ii0)>2
%             P=ii0-min_ind;
%             [~,I]=sort(P,'ascend'); %小さい順に並び替え
% %             ii=[ii(I(1)) ii(I(2))];
%             zi0=[zi0(I(1)) zi0(I(2))];
% %             yi=[yi(I(1)) yi(I(2))];
%         end
%         zi(:,i)=zi0;
%         ri(:,i)=polyval(v(:,i),[zi0(1) zi0(2)]);
%         clear zi0 yi0 ii0


        %% ver1:電流シート幅の定義（Jt最大の点）
%         sheetJt(i)=min(Jt(:,i)).*0.5;　%(5-1:poly)垂線上のJtとmin(jt)*0.5の交点による半値全幅FWHM
%         sheetJt(i)=0;　%(5-2:fw)垂線上のJtとJt=0の交点による全値幅FW
%         sheetJt(i)=Jtmin(:,:,i).*0.5;　%(5-3:fwhm03)全rz領域での最小のJtに対する半値全幅
        sheetJt(i)=-1e+5; %(5-4:fw04)閾値Jt=-1e+5を設定した全値幅FW

        sheetJtarray=ones(nq,1).*sheetJt(i);

        [zi0,yi0,ii0]=polyxpoly(grid2D.zq(1,1:nq),Jt(1:nq,i),grid2D.zq(1,1:nq),sheetJtarray);
        if size(ii0,1)>2
            P=ii0(:,2)-min_ind.*ones(size(ii0,1),1);
            P=abs(P);
            [~,I]=sort(P,'ascend'); %小さい順に並び替え
%             ii=[ii(I(1)) ii(I(2))];
            zi0=[zi0(I(1)) zi0(I(2))];
%             yi=[yi(I(1)) yi(I(2))];
        end
        zi(1:min(size(ii0,1),2),i)=zi0;
        ri(1:min(size(ii0,1),2),i)=polyval(v(:,i),zi0);
        if size(ii0,1)<2
            zi(size(ii0,1)+1:2,i)=NaN;
            ri(size(ii0,1)+1:2,i)=NaN;
        end
        
        %シート幅計算
        if size(ii0,1)==2
            wide(i)=norm([zi(2,i)-zi(1,i) ri(2,i)-ri(1,i)]);
        else
            wide(i)=NaN;
        end
        clear zi0 yi0 ii0

        %% ver2:電流シート幅の定義 :Br(=Bp)の極から求める場合
    max_log=islocalmax(Bp(:,i));
    max_inf=find(max_log);
    min_log=islocalmin(Bp(:,i));
    min_inf=find(min_log);
    if numel(max_inf)==0
        zi(2,i)=NaN;
        zi(1,i)=NaN;
        ri(2,i)=NaN;
        ri(1,i)=NaN;
        wide(i)=NaN;
    elseif numel(min_inf)==0
        zi(2,i)=NaN;
        zi(1,i)=NaN;
        ri(2,i)=NaN;
        ri(1,i)=NaN;
        wide(i)=NaN;
    else
        zi(2,i)=grid2D.zq(1,max_inf);
        zi(1,i)=grid2D.zq(1,min_inf);
        ri(2,i)=vertical(max_inf);
        ri(1,i)=vertical(min_inf);
        wide(i)=norm([zi(2,i)-zi(1,i) ri(2,i)-ri(1,i)]);
    end

    end

    %磁気面プロット上に表示して確認
    subplot(row,column,i)
    contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),8,'w-','LineWidth',0.1);
    colormap(jet)
    caxis([-2.5*1e+6,1.6*1e+6])
    axis image
    axis tight manual
    hold on
    contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black','LineWidth',0.1)
%     colorbar('Location','eastoutside')
    plot(diff,grid2D.rq(r_s:r_n,1),'b-','LineWidth',1)
    plot(grid2D.zq(1,:),vertical,'b-','LineWidth',1)
    plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"bx",'MarkerSize',8)%Xpoint
%     plot(grid2D.zq(1,z_xind(:,i)),grid2D.rq(r_xind(:,i),1),"bx",'MarkerSize',8)%Xpoint
    plot(zi(1,i),ri(1,i),'b|','MarkerSize',12)
    plot(zi(2,i),ri(2,i),'b|','MarkerSize',12)
    plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1),'b-','LineWidth',0.5)
    hold off
    alpha(0.3)
%     xlabel('z [m]')
%     ylabel('r [m]')
    title(string(tstart+i-1)+' us')
end
sgtitle(strcat('IDX=',num2str(IDX)))
% saveas(gcf,strcat(pathname.save,'\',num2str(IDX),'wide_Br.png'))
% close

%垂線上のJtのplot
f2=figure;
f2.WindowState='maximized';
for i=1:tsize
     subplot(row,column,i)
     plot(grid2D.zq(1,:),Bp(:,i))%plot(grid2D.zq(1,:),Jt(:,i))
     xline(grid2D.zq(1,z_xind(:,i)),'r--')
%      yline(sheetJt(i),'b--')%yline(min(Jt(:,i))*0.5,'b--')
     xline(zi(1,i),'b-')%xline(grid2D.zq(1,left_ind),'b-')
     xline(zi(2,i),'b-')%xline(grid2D.zq(1,right_ind),'b-')
     yline(0,'k-')
     xlabel('z [m]')
     ylabel('Br [Wb]')%ylabel('Jt [A/m^{2}]')
     xlim([grid2D.zq(1,1) grid2D.zq(1,n)])
%      ylim([-1.6e+6 0.1e+6])%(a)
%      ylim([-1e+6 0.1e+6])%(b)
%      ylim([-0.5e+6 0.1e+6])%(c)
     title(string(tstart+i-1)+' us')
end
sgtitle(strcat('IDX=',num2str(IDX)))
% saveas(gcf,strcat(pathname.save,'\',num2str(IDX),'Br.png'))
% close
     
end



        






