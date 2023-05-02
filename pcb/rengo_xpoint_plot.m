%% 
% X点でのpsi_commonや電流密度、電場を時系列で求めるコード
% 
% plot開始点は合体率が0を超える時間から
% 
% （コード重いので3こずつくらいで出力した方が安心）

clearvars
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス
pathname.fourier=getenv('fourier_path');%md0までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath
pathname.save=getenv('save_path'); %保存先

DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';
T=getTS6log(DOCID);% ログのテーブルを取得

IDXlist=[2911:2913 2925 2926 2927 2931 2933 2947:2950 2942 2943 2946];

start=465;
frame=45;
for i=1 %numel(IDXlist)
IDX=IDXlist(i);
tstart=472;
tend=482;
xdata(T, IDX, pathname,start,frame,tstart,tend);
end

function xdata(T, IDX, pathname,start,frame,tstart,tend)
date=T.date(IDX);
shot=T.shot(IDX);
TF_shot=T.TFoffset(IDX);
d_tacq=T.d_tacq(IDX);
d_tacqTF=T.TFdtacq(IDX);
offset_TF=isfinite(TF_shot);
Cr=T.Crowbar_us_(IDX);
if isnan(Cr)
    Cr=num2str('-');
end
CB=T.CB1_kV_(IDX);
%start=T.Period_StartTime_(IDX);
if isnan(T.EF_A_(IDX))
    i_EF=150;
else  
    i_EF=T.EF_A_(IDX);
end
n=50;

%CAO-probe(Et,Jt,merge ratio)
[B_z,r_probe,z_probe,ch_dist,data,data_raw,shot_num]= get_B_z(date,TF_shot,shot,offset_TF,i_EF,pathname.ts3u);
B_z = B_z([2,3,4,6,7,8],2:end,:);
data = data([2,3,4,6,7,8],2:end,:);
z_probe = z_probe(2:end);
ch_dist = ch_dist([2,3,4,6,7,8],2:end);
r_probe = r_probe([2,3,4,6,7,8]);

mesh_z=100;
mesh_r=100;

z_space = linspace(z_probe(1),z_probe(end),mesh_z);
r_space = linspace(r_probe(1),r_probe(end),mesh_r);
[psi_mesh_z,psi_mesh_r] = meshgrid(z_space,r_space);

psi_store = zeros(length(r_space),length(z_space),frame);
Jt_store = psi_store;
Et_store = psi_store;
for i=1:frame
    psi  = get_psi(B_z,r_probe,i+start);
    psi_store(:,:,i)= griddata(z_probe,r_probe,psi,psi_mesh_z,psi_mesh_r,'cubic');
    [j_t,z_space_jt,r_space_jt] = jt(B_z,z_probe,r_probe,i+start);
    Jt_store(:,:,i) = griddata(z_space_jt,r_space_jt,j_t,psi_mesh_z,psi_mesh_r,'v4');
    E_t = Et(B_z,r_probe,i+start);
    Et_store(:,:,i) = griddata(z_probe,r_probe,E_t,psi_mesh_z,psi_mesh_r,'v4');
end

[max_psi,max_psi_r]=max(psi_store,[],1); %各時間、各列(z)ごとのpsiの最大値
max_psi=squeeze(max_psi);
psi_pr=zeros(3,frame);
xJt=zeros(1,frame);
xEt=xJt;
%plot(max_psi)
for i=1:frame
    r_ind=max_psi_r(:,:,i);
    max_psi_ind=find(islocalmax(smooth(max_psi(:,i)),'MaxNumExtrema', 2));
    if numel(max_psi(min(max_psi_ind),i))==0
    psi_pr(1,i)=NaN;
    else
    psi_pr(1,i)=max_psi(min(max_psi_ind),i);
    end
    if numel(max_psi(max(max_psi_ind),i))==0
    psi_pr(2,i)=NaN;
    else
    psi_pr(2,i)=max_psi(max(max_psi_ind),i);
    end
    if numel(find(islocalmin(smooth(max_psi(:,i)),'MaxNumExtrema', 1)))==0
        psi_pr(3,i)=NaN;
        xJt(1,i)=NaN;
        xEt(1,i)=NaN;
    else
        min_psi_ind=islocalmin(smooth(max_psi(:,i)),'MaxNumExtrema', 1);
        xr=r_ind(min_psi_ind);
        if xr==1 || xr==mesh_r %r両端の場合は検知しない
            psi_pr(3,i)=NaN;
            xJt(1,i)=NaN;
            xEt(1,i)=NaN;
        else
        psi_pr(3,i)=max_psi(min_psi_ind,i);
        xJt(1,i)=Jt_store(xr,min_psi_ind,i);
        xEt(1,i)=Et_store(xr,min_psi_ind,i);
        end
    end
    if max_psi(1,i)==max(max_psi(1:end/2,i))
        psi_pr(1,i)=NaN; %max_psi(1,i);
    end    
    if max_psi(end,i)==max(max_psi(end/2:end,i))
        psi_pr(2,i)=NaN; %max_psi(end,i);
    end
end
fitrate=psi_pr(3,:)./min(psi_pr(1:2,:),[],1);
merge_ind=find(fitrate>=0,1); %fitrateが初めて0以上となる点のindex

time=[1:frame]+start;

%pcb probe(slope)
trange = tstart:tend+1;
tsize = tend-tstart+1;
%[date, ~, ~, ~, i_EF, ~, ~, d_tacq, d_tacqTF, n] = getinput(T,IDX);%Tのテーブルから入力のリストを出力
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
for i=1:tsize
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
        %diff=polyval(p(:,i),grid2D.rq(r_s:r_n,1));
        if z_xind==1 || z_xind==n
            p(:,i)=NaN;
        end
    else
        p(:,i)=NaN;
    end   
end

%X点での各データのplot
figure('Position',[0,0,600,400])
subplot(3,1,1)
plot(time,fitrate,'-o','Color','k','LineWidth',1)
xline(472,'b-',"LineWidth",1)
xline(475,'r-',"LineWidth",1)
%yline(1,'k')
ylim([0 1])
ylabel('Merging ratio','FontSize',13)
%xlim([time(merge_ind) 486])
%xlabel('Time [us]')
xlim([time(merge_ind) time(merge_ind)+20])
ha1 = gca;
ha1.LineWidth = 1;
ha1.FontSize=13;
ha1.XTickLabel = cell(size(ha1.XTickLabel));

subplot(3,1,2)
yyaxis left
plot(time,-1.*xJt.*1e-6,'-o','Color','b','LineWidth',1)
%yline(0,'k-')
xline(472,'b-',"LineWidth",1)
ylabel('Jt at Xpoint [MA/m^{2}]','FontSize',13,'Color','b')
xlabel('Time [us]','FontSize',13)
%ylim([-10e+5 3e+5])
ylim([0 1.1])
xlim([time(merge_ind) time(merge_ind)+20])
yyaxis right
plot(time,-1.*xEt.*1e-3,'-o','Color','r','LineWidth',1)
xline(475,'r-',"LineWidth",1)
%yline(0,'k-')
ylabel('Et at Xpoint [kV/m]','FontSize',13,'Color','r')
xlabel('Time [us]','FontSize',13)
%ylim([-3.3e-4 0.5e-4])
ylim([0 0.35])
xlim([time(merge_ind) time(merge_ind)+20])
%set(gca,'FontSize',13);
ha2 = gca;
ha2.LineWidth = 1;
ha2.FontSize=13;

subplot(3,1,3)
plot(tstart:tend,p(1,:).*180./pi,'ko-','LineWidth',1)
xlim([time(merge_ind) time(merge_ind)+20])
ylim([-0.55 0.1].*180./pi)
xlabel('Time [us]')
ylabel('sheet angle')
ha3 = gca;
ha3.LineWidth = 1;
ha3.FontSize=13;

sgtitle(strcat('IDX',num2str(IDX),', Cr',num2str(Cr),' us, CB',num2str(CB),' kV'),'FontSize',12)
%saveas(gcf,strcat(pathname.save,'\IDX',num2str(IDX),'_xdata.png'))


end
