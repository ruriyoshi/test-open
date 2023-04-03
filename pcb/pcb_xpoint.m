%X点での電場、Jtを求める
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス
pathname.fourier=getenv('fourier_path');%md0までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath
pathname.save=getenv('save_path');

DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';
T=getTS6log(DOCID);% ログのテーブルを取得

IDXlist=[2911:2913 2925 2926 2927 2931 2933];
tstartlist=[472 472 472 477 477 477 482 482];
tendlist=[482 482 482 487 487 487 500 500];

i=8;
IDX=IDXlist(i);
tstart=477;
tend=496;
[Jt,Et,eta]=cal_xpoint(T,IDX, pathname, tstart, tend);

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

function [Jt,Et,eta]=cal_xpoint(T,IDX, pathname, tstart, tend)
trange = tstart:tend+1;
tsize = tend-tstart+1;

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
Jt=zeros(tsize,1);
Et=zeros(tsize,1);
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
        Jt(i)=data2D.Jt(r_xind,z_xind,i);
        Et(i)=data2D.Et(r_xind,z_xind,i);
        if z_xind==1 || z_xind==n
            Jt(i)=NaN;
            Et(i)=NaN;
        end
    else
        Jt(i)=NaN;
        Et(i)=NaN;
    end   
end
eta=Et./Jt;

f2=figure;
f2.WindowState='maximized';
subplot(3,1,1)
plot(tstart:tend,Jt.*1e-6,'b+-')
xlim([tstart tend])
ylabel('Jt [MA/m^{2}]')
subplot(3,1,2)
plot(tstart:tend,Et.*1e-3,'r+-')
xlim([tstart tend])
ylabel('Et [kV/m]')
subplot(3,1,3)
plot(tstart:tend,eta.*1e+3,'r+-')
xlim([tstart tend])
ylabel('η [mΩ・m]')
xlabel('Time [us]')
ha3 = gca;
ha3.LineWidth = 1;
ha3.FontSize=13;
sgtitle(strcat('IDX=',num2str(IDX)))
% saveas(gcf,strcat(pathname.save,'\',num2str(IDX),'xpoint.png'))

end