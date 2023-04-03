%電流シートの傾きを求める（完成版）
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス
pathname.fourier=getenv('fourier_path');%md0までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath
pathname.save=getenv('save_path');

DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';
T=getTS6log(DOCID);% ログのテーブルを取得

%IDXlist=[2911:2913 2925 2926 2927 2931 2933 2947:2950 2942 2943 2946];
% IDXlist=[2911 2926 2931 2948 2942];
% tstartlist=[472 477 482 477 484];
% tendlist=[482 488 503 487 499];

IDXlist=[2911:2913 2925 2926 2927 2931 2933];
tstartlist=[472 472 472 477 477 477 482 482];
tendlist=[482 482 482 487 487 487 500 500];

figure
aslope=zeros(3,11);
for i=1:3
    IDX=IDXlist(i);
    tstart=tstartlist(i);
    tend=tendlist(i);
    p=cal_diff(T, IDX, pathname, tstart, tend);
    plot(tstart:tend,p(1,:).*180./pi,'ko-','LineWidth',1)
    aslope(i,:)=p(1,:).*180./pi;
    hold on
end
bslope=zeros(3,11);
for i=4:6
    IDX=IDXlist(i);
    tstart=tstartlist(i);
    tend=tendlist(i);
    p=cal_diff(T, IDX, pathname, tstart, tend);
    plot(tstart:tend,p(1,:).*180./pi,'ro-','LineWidth',1)
    bslope(i-3,:)=p(1,:).*180./pi;
    hold on
end
cslope=zeros(2,19);
for i=7:8
    IDX=IDXlist(i);
    tstart=tstartlist(i);
    tend=tendlist(i);
    p=cal_diff(T, IDX, pathname, tstart, tend);
    plot(tstart:tend,p(1,:).*180./pi,'bo-','LineWidth',1)
    cslope(i-6,:)=p(1,:).*180./pi;
    hold on
end
yline(0,'k-')
hold off
xlim([470 502])
xlabel('Time [us]')
ylabel('maximum sheet angle [degree]')
ylim([-0.55 0.1].*180./pi)
%legend('IDX2911','IDX2926','IDX2931','Location','southeast')
ha3 = gca;
ha3.LineWidth = 1;
ha3.FontSize=13;

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

function p=cal_diff(T, IDX, pathname, tstart, tend)
trange = tstart:tend+1;
tsize = tend-tstart+1;
row = 2;
column = ceil(tsize/row); 

[date, ~, ~, ~, i_EF, ~, ~, d_tacq, d_tacqTF, n] = getinput(T,IDX);%Tのテーブルから入力のリストを出力
[grid2D, data2D] = pcbdata(date, d_tacq, d_tacqTF, trange, [], n,i_EF);

if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end

%%%midplaneとかO点、X点を探す
[psimid,mid]=min(data2D.psi,[],2);
[opoint,~]=islocalmin(psimid,1);
[xpoint,~]=islocalmax(psimid,1);
% [xp_psi,maxxp]=max(squeeze(psimid),[],1);
r_index=1:n;

%X点での勾配
p=zeros(2,tsize);
% figure('Position', [0 0 1500 1500],'visible','on');
for i=1:tsize
%     r_xind=ceil(n/2);
%     z_xind=ceil(n/2);
    if sum(xpoint(:,:,i))>0   
        r_xind=r_index(xpoint(:,:,i));
        z_xind=mid(xpoint(:,:,i),:,i);
        if numel(r_xind)>1
            [~,I]=max(psimid(r_xind));
            r_xind=r_xind(I);
            z_xind=z_xind(I);
        end
        
        r_s=max(2,r_xind-5);
        r_n=min(n,r_xind+5);
        p(:,i)=polyfit(grid2D.rq(r_s:r_n,1),grid2D.zq(1,squeeze(mid(r_s:r_n,:,i))),1);
%         diff=polyval(p(:,i),grid2D.rq(r_s:r_n,1));
%         subplot(row,column,i)
%         contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20)
%         axis image
%         axis tight manual
%         hold on
%         plot(diff,grid2D.rq(r_s:r_n,1),'r-','LineWidth',1)
%         plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"rx")
%         plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1),'b-')
%         hold off
%         title(string(tstart+i-1)+' us')

        if z_xind==1 || z_xind==n
            p(:,i)=NaN;
        end
    else
        p(:,i)=NaN;
    end   
end
% sgtitle(strcat('IDX=',num2str(IDX)))

% figure
% plot(tstart:tend,p(1,:),'ko-')
% xlim([tstart tend])
% xlabel('Time [us]')
% ylabel('current sheet slope')
% ylim([-0.55 0.1])
% ha3 = gca;
% ha3.LineWidth = 1;
% ha3.FontSize=13;
% Cr=T.Crowbar_us_(IDX);
% if isnan(Cr)
%     Cr=num2str('-');
% end
% CB=T.CB1_kV_(IDX);
% sgtitle(strcat('IDX',num2str(IDX),', Cr',num2str(Cr),' us, CB',num2str(CB),' kV'),'FontSize',12)
end