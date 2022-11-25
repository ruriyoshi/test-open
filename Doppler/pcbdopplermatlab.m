clear all
close all
clear
%%%%%%%%%%%%%%%%%%%%%%%%
%  288chDopplerとpcbデータを重ねてプロットして保存するコード
%　288chDopplerはIDLで作ったsavのデータをあらかじめI:\makimitsu\yyddmmに保存してあるものを読み込む
%
%%%%%%%%%%%%%%%%%%%%%%%%

yourname = 'C:\Users\Moe Akimitsu\';
f = fullfile(yourname,'Documents','GitHub','test-open');
addpath(genpath(f));
%addpath(fullfile(yourname,'Documents','GitHub','test-open','pcb'));
%f = fullfile(yourname,'Documents','GitHub','SXR_test');
%addpath(f);
%%%適宜変更


%DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';
DOCID='1spjcYNNser4x_oR03uNuyLKBG-eXmbyul9-I5UDkiuI';
T=getTS6log(DOCID);% ログのテーブルを取得
 
%shotlist=[2694:2754];
%shotlist=[2785:2817];
shotlist=[50];
subT=T(shotlist,:);
IDXlist=shotlist(isfinite(subT.DopplerDelay)&isfinite(subT.d_tacq));
%IDX=IDXlist(1,88);
for IDX=IDXlist%(1,49:end)  %
    close all
date=T.date(IDX)
shot=T.shot(IDX)
d_tacq=T.d_tacq(IDX);
d_tacqTF=T.TFdtacq(IDX);
if isnan(T.EF_A_(IDX))%%NaNでないことを確認（ログが空白だとNaNになる）
    i_EF=150;
else  %NaNなら150をとりあえず代入、記入されているときはその値を使う
    i_EF=T.EF_A_(IDX);
end
trange=460:490;
t=T.DopplerDelay(IDX);
%%%pcbの部分
 n=30;
 [grid2D, data2D] = pcbdata(date, d_tacq,d_tacqTF,trange, [], n,i_EF);
% if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
%     continue
% end
%     maxrange=max(abs(data2D.Jt),[],'all');
% [psimid,mid]=min(data2D.psi,[],2);
% [opoint,p]=islocalmin(psimid,1);
% [xpoint,~]=islocalmax(psimid,1);
% [xp_psi,maxxp]=max(squeeze(psimid),[],1);
% onum=squeeze(sum(opoint,1));
% trange(onum~=0)
    %maxrange=2e6;

%%磁気面時間発展プロット
% figure
% start=8; %460+?
% for m=1:10 
%     i=start+m;
%     t=trange(i);
%     subplot(2,5,m)
%     contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Jt(:,:,i),10,'LineStyle','none')
%     colormap(jet)
%     axis image
%     axis tight manual
%     %     xlim([-0.02 0.02])
%     %     ylim([0.12 0.27])
%     caxis([-maxrange,maxrange])
%     colorbar('Location','eastoutside')
%     %zlim([-1 1])
%     %colormap(bone)
%     hold on
%     plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
%     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),50,'black')
%     plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"ro")
%     plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"rx")
%     hold off
%     title(string(t)+'us')
%     xlabel('z')
%     ylabel('r')
% end

offsetsmile=0;
%%ファイルへのパスを作る
%それぞれのPCから共有フォルダまでのパスはそれぞれ異なるので各自で設定
%pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス
pathname.fourier='I:';%md0までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath
filenameasc= dir(strcat(fullfile(pathname.fourier,'Doppler','288CH',strcat('20',num2str(date))),'\shot',num2str(shot),'_*.asc'));
if numel(filenameasc)==1%共有フォルダ以下から目的ショットのファイルを探す
filenameasc=fullfile(filenameasc.folder, filenameasc.name);
tmax=300;
EM_max=1.e8;
lambda0=480.6;%656.3%471.3%656.3%480.6%587.562%486.133%656.3%486.133%468.57%486.133%nm
mass=39.95;%4.%12.01%39.95%1.
resolution=-1.420387e-12*lambda0^3 - 2.156031e-09*lambda0^2 + 1.250038e-06*lambda0 + 3.830769e-03;%0.0037714% -0.000000000001420387*lambda0^3 - 0.000000002156031*lambda0^2 + 0.000001250038*lambda0 + 0.003830769

z=readmatrix("z.txt")*1e-3;
%[file,path] = uigetfile('*.asc')
%filename=fullfile(path,file)
[d,delimiterOut]=importdata(filenameasc,'	');d=d(:,2:end)';
%[file,path] = uigetfile('*.asc')
%filenamebg=fullfile(path,file)

%calib=readmatrix("Ar.calibration.0604.txt");
if date==211223 || date==211224
filenamebg='I:\Doppler\288CH\20211224\bg.asc';
calib=readmatrix("Ar_calibration_20211223-24.txt");
elseif date==211217
filenamebg='I:\Doppler\288CH\20211217\bg.asc';
calib=readmatrix("Ar_calibration_20211217.txt");
elseif date==211219
filenamebg='I:\Doppler\288CH\20211219\bg.asc';
calib=readmatrix("H_calibration_20211219.txt");
lambda0=486.133;%656.3%471.3%656.3%480.6%587.562%486.133%656.3%486.133%468.57%486.133%nm
mass=1;%4.%12.01%39.95%1.
resolution=-1.420387e-12*lambda0^3 - 2.156031e-09*lambda0^2 + 1.250038e-06*lambda0 + 3.830769e-03;%0.0037714% -0.000000000001420387*lambda0^3 - 0.000000002156031*lambda0^2 + 0.000001250038*lambda0 + 0.003830769
elseif date==190608
filenamebg='I:\Doppler\288CH\20190604\bg.asc';
calib=readmatrix("Ar.calibration.0604.txt");  
    if shot<=39
        offsetsmile=50;
    else
        offsetsmile=-5;
    end    
end

bg=importdata(filenamebg);bg=bg(:,2:end)';d=d-bg;

CH=squeeze(calib(:,1));
Center=squeeze(calib(:,2))-2;
Smile=squeeze(calib(:,3))+offsetsmile;
relative=squeeze(calib(:,5));
instru=squeeze(calib(:,6));


separation=find(mod(CH-1,16)==1)-1 ; separation(4)=separation(4)+1;
%なんかあわない

p=readmatrix("r.txt")*1.e-3;p=p(:,4)
zall=repmat(z,1,16)';
pall=repmat([p,flip(p)],1,9);
edge=0.3;
Ti_instru=1.69e8.*mass.*(2.*resolution.*instru.*sqrt(2.*log(2.))/lambda0).^2;

lambda=zeros([1024 numel(CH)]);
x=(1:1024);
for i=1:numel(CH) 
lambda(:,i)=(x-Smile(i))*resolution+lambda0;%+0.13;
end
spectra=zeros([1024 numel(CH)]);
figure('Visible','on')
subplot('Position',[0.1,0.4,0.85,0.55])
interval=7;
contourf(x,x,d,30,'LineStyle','none')
%caxis([min(d,[],'all') (max(d,[],'all')-min(d,[],'all'))*0.1+min(d,[],'all')])
caxis([min(d,[],'all') max(d,[],'all')])
ylim([min(Center)-interval-5 max(Center)+interval+5])
axis tight
hold on
for i=1:numel(CH)
plot(Smile(i)+[-40,40],Center(i)*[1,1],'r')
plot(Smile(i)*[1.,1.],Center(i)+[-7,7],'g')
plot(Smile(i)+[-40,40],Center(i)*[1,1]+interval,':y')
plot(Smile(i)+[-40,40],Center(i)*[1,1]-interval,':b')
end
hold off

 xr=[410,540]
 lambdaA=(x-(xr(1)+xr(2))/2.).*resolution+lambda0 ;
 lambdaB=lambdaA 
 lambdaA=lambdaA(xr(1):xr(2));
ybin=sum(d,1);
subplot('Position',[0.1,0.05,0.85,0.35])
plot(x,ybin,'*')
axis tight
hold on
%gaussEqn='a0*exp(-((x-a1)/a2)^2)+a3+a4*x';
f=fit(x(xr(1):xr(2))',ybin(xr(1):xr(2))','gauss1')
plot(f,x(xr(1):xr(2)),ybin(xr(1):xr(2)))

% A=[coeff[0],200,coeff[2],coeff[0],350,coeff[2],coeff[0],500,coeff[2],coeff[0],650,coeff[2],coeff[0],800,coeff[2],coeff[3],0,0]
% fita=A*0.+1.
% %stop
% %**************************
% %*Line-integrated analysis*
% %**************************
% 
 passive.Ti=CH*0.;
 passive.Timax=CH*0.;
 passive.Timin=CH*0.;
 passive.Em=CH*0.; 

%figure('Visible','on')
for i=1:numel(CH)
%     if ismember(i,separation)
%         figure('Name',num2str(round(CH(i)/16+2)),'Visible','on')
%     end
%    subplot(4,4,mod(CH(i)-1,16)+1)
    spectra(:,i)=sum(d(Center(i)-interval:Center(i)+interval,:),1)*relative(i)    ;    %各チャンネルのガウシアン信号を積分
    input=squeeze(spectra(:,i))-median(spectra(:,i)) ;
    c=fit(x(1,round(Smile(i))-75:round(Smile(i))+75)',input(round(Smile(i))-75:round(Smile(i))+75,1),'gauss1') ;  %passive.Ti(i)=1.69e8*mass*(2.*resolution*(coeff[2])*sqrt(2.*alog(2.))/lambda0)^2-Ti_instru(i)
    conf=confint(c);
    passive.Ti(i)=1.69e8*mass*(2.*resolution*c.c1*sqrt(log(2.))/lambda0)^2-Ti_instru(i);
    passive.Timax(i)=1.69e8*mass*(2.*resolution*conf(2,3)*sqrt(log(2.))/lambda0)^2-Ti_instru(i);
    passive.Timin(i)=1.69e8*mass*(2.*resolution*conf(1,3)*sqrt(log(2.))/lambda0)^2-Ti_instru(i);
%    plot(c,x(1,round(Smile(i))-75:round(Smile(i))+75),input(round(Smile(i))-75:round(Smile(i))+75,1))
    passive.Em(i)=resolution*sum(input(round(Smile(i))-75:round(Smile(i))+75));
%    spectra(:,i)=input;
    lambda(:,i)=(x-c.b1)*resolution+lambda0;
end
clear conf c input 


figure('Visible','on')
plot(CH,passive.Ti,'ro')
hold on
plot([CH CH]',[passive.Timin passive.Timax]','b-_')
xline(CH(separation),':')
%ar1=area(CH,[passive.Timax, passive.Timin-passive.Timax])
%set(ar1(1),'FaceColor','None','LineStyle','None')
%set(ar1(2),'FaceColor',[0 0.2 1],'FaceAlpha',0.2,'LineStyle','None')
% errplot,CH,passive.Timin,


figure('Visible','on')
 Ti2D=zeros([numel(p),numel(z)]);
 Em2D=zeros([numel(p),numel(z)]);
 for i=1:size(zall,2)
 Ti2D(:,i)=pchip(pall(CH(zall(CH)==zall(1,i))),passive.Ti(zall(CH)==zall(1,i)),sort(pall(:,i)));
 Em2D(:,i)=pchip(pall(CH(zall(CH)==zall(1,i))),passive.Em(zall(CH)==zall(1,i)),sort(pall(:,i)));
 plot(pall(CH(zall(CH)==zall(1,i))),passive.Em(zall(CH)==zall(1,i)),'o')
 hold on
 plot(sort(pall(:,i)),Em2D(:,i))
 end
 % passiveGrid.zq=zall;
% passiveGrid.rq=pall;

 figure('Visible', 'on')
 contourf(z,p,(Ti2D)/max(Ti2D,[],'all'),60,'LineStyle','none')
 figure('Visible', 'on')
 contourf(z,p,(Em2D)/max(Em2D,[],'all'),60,'LineStyle','none')

% %******************
% %*Abel-ininversion*
% %******************

num=100;
%yy=[1:num]/num*(edge-min(p))+min(p);
yy=[1:num]/num*(edge-0.06)+0.06;
spectra_interp=zeros([numel(yy),numel(lambdaA),numel(z)]);
%grid2D=meshgrid(lambdaA,yy);
for i=1:numel(z) 
    [Yin,I]= sort(pall(CH(zall(CH)==zall(1,i))));
    Xin=lambda(:,zall(CH)==zall(1,i));
    Zin=spectra(:,zall(CH)==zall(1,i));
    Yin=[0.06, Yin', edge];
    Xin=[lambdaB',Xin(:,I),lambdaB'];
    Zin=[zeros(numel(lambdaB),1),Zin(:,I),zeros(numel(lambdaB),1)];
    result=Trigrid_interpor_for_r_lambda(Zin,Xin,Yin,lambdaA,yy);
spectra_interp(:,:,i)=result.Z;
end
% for j=1:numel(z)
%  for i=1:num  
%    if mod(i, 36)==1 
%        figure('Visible','on')
%    end
%    subplot(6,6,mod(i-1, 36)+1)
%    f=fit(lambdaA',squeeze(spectra_interp(i,:,j))','gauss1');
%    plot(f,lambdaA',squeeze(spectra_interp(i,:,j)'))
%    legend('off')
%  end 
% end

 
%  spectra_interp=smoothdata(spectra_interp,1,'movmean',num/16);
% spectra_interp=smoothdata(spectra_interp,2,'movmean',5);
% spectra_interp=smoothdata(spectra_interp,3,'movmean',1);
% dy=yy(2)-yy(1);
reconst_int=spectra_interp*0;
[~,dndy,~]=gradient(spectra_interp,1,yy,1);
dy=diff(yy);
for i=1:numel(z)
 for j=1:num 
  for k=j+1:num 
     reconst_int(j,:,i)=reconst_int(j,:,i)-1./pi.*dndy(k,:,i)./sqrt(yy(k).^2-yy(j).^2).*dy(k-1);
  end
 end
end
% reconst_int=smooth(reconst_int,[5,num/16,1])
% 
 emission=zeros([num,numel(z)]);
 Ti.Ti2D=zeros([num,numel(z)]);
 Ti.max=zeros([num,numel(z)]);
 Ti.min=zeros([num,numel(z)]);
 Ti_instru2=Ti_instru*0;
 for i=1:numel(z)
     Ti_instru2(zall(CH)==zall(1,i))=sum(Ti_instru(zall(CH)==zall(1,i)))/sum(zall(CH)==zall(1,i));
 end 
 for j=1:numel(z)
% figure('Visible','on')
    for i=1:num  

%    if mod(i, 36)==1 
%        figure('Visible','on')
%    end
%  subplot(ceil(sqrt(num)),ceil(num/sqrt(num)),i)
    input=squeeze(reconst_int(i,:,j));%input(input<=0)=-input(input<=0);
   if sum(input)<=0
       checker=false;
   else     
    c=fit(lambdaA',input','gauss1','StartPoint',[max(input),480.6,0.002]);   
    conf=confint(c);
    Ti.Ti2D(i,j)=1.69e8*mass*(2.*c.c1*sqrt(log(2.))/lambda0)^2-Ti_instru2(j);
    Ti.max(i,j)=1.69e8*mass*(2.*conf(2,3)*sqrt(log(2.))/lambda0)^2-Ti_instru2(j);
    Ti.min(i,j)=1.69e8*mass*(2.*conf(1,3)*sqrt(log(2.))/lambda0)^2-Ti_instru2(j);
    emission(i,j)=resolution*sum(input);
   checker=(abs(c.b1-lambda0) < 0.3) & (c.a1 > 0 )&(emission(i,j) >= 100);
   % checker=(abs(c.b1-lambda0) < 0.3) & (abs(Ti.max(i,j)-Ti.min(i,j)) < Ti.Ti2D(i,j))%*float(abs(Ti.max(i,j)-Ti.min(i,j)) lt Ti.Ti2D(i,j)+Ti_instru2[j])%*(emission_local(i,j) gt EM_max*0.1);
   end
 if checker==false
     Ti.Ti2D(i,j)=NaN;
     Ti.max(i,j)=NaN;
     Ti.min(i,j)=NaN;
     emission(i,j)=NaN;
   %  plot(lambdaA',input',':')
 else 
   % plot(c,lambdaA',input')
 end
 end 
 end  
doppler=struct('z',z,'yy',yy,'emission',emission,'ti_2d',Ti.Ti2D);
    i=t-459;
%    restore_idl( fullfile(filepath.D288.folder,filepath.D288.name),'lowercase','create'); %.savファイルを読み込む
    
    f = figure;
    % f.WindowState = 'maximized';
    f.Units = 'normalized';
    f.Position = [0.1,0.2,0.8,0.6];
    pos1 = [0.07,0.2,0.35,0.6];
    pos2 = [0.58,0.2,0.35,0.6];

    subplot('Position',pos1);
    contourf(doppler.z,doppler.yy,doppler.emission,30,'LineStyle','none')
    colormap(jet)
    axis image
    axis tight manual
    c=colorbar('Location','eastoutside');
    c.Label.String = 'Emission [a.u.]';
    hold on
    plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
    contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),30,'black')
    plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"ro")
    plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"rx")
    hold off
    title(string(t)+'us,emiision')
    xlabel('z')
    ylabel('r')
     xlim([-0.05 0.05])
    ylim([0.07 0.25])
    caxis([-3e5,3e5])

    subplot('Position',pos2);
    contourf(doppler.z,doppler.yy,doppler.ti_2d,[0:2:150],'LineStyle','none')
    colormap(jet)
    axis image
    axis tight manual
    colorbar('Location','eastoutside')
    hold on
    plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
    contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),30,'black')
    plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"ro")
    plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"rx")
    hold off
    %caxis([-100,100])
    title(string(t)+'us,ti')
    xlim([-0.05 0.05])
   ylim([0.07 0.25])
    xlabel('z')
    ylabel('r')
    cd 'C:\Users\Moe Akimitsu\Desktop\dronkaiseki3\'
    %filename = strcat('J:\makimitsu\',num2str(date),'\Doppler_',num2str(date),num2str(shot,'%03i'),'_',num2str(t),'us');
    %filename = strcat('C:\Users\Moe Akimitsu\Desktop\dronkaiseki\',num2str(date),'\Doppler_',num2str(date),num2str(shot,'%03i'),'_',num2str(t),'us');
                  filename=num2str(IDX);
    saveas(gcf,strcat(filename,'.png'))
    close
end


% %[B_z,r_probe,z_probe,ch_dist,B_z_return,data_return,shot_num] 
% [B_z,r_probe,z_probe,ch_dist,data,data_raw,shot_num]= get_B_z(date,TF_shot,shot,offset_TF,T.EF_A_(IDX),pathname.ts3u);
% B_z = B_z([2,3,4,6,7,8],2:end,:);
% data = data([2,3,4,6,7,8],2:end,:);
% z_probe = z_probe(2:end);
% ch_dist = ch_dist([2,3,4,6,7,8],2:end);
% r_probe = r_probe([2,3,4,6,7,8]);



%plot_B_z_in_time(B_z,ch_dist,350,600);
%plot_psi_SXR_multi(B_z,r_probe,z_probe,date,shot,layer,area,start,exposure,SXRfilename)
%plot_psi_SXR_multi(B_z,r_probe,z_probe,date,shot,true,true,T.Period_StartTime_(IDX),2,filepath.SXR)
end






 function result=Trigrid_interpor_for_r_lambda(z,x,y,xout,yout)

%Nterms=N_elements(x)*N_elements(y)
Nterms=numel(x);
Xtemp=x' ;Xtemp= Xtemp(:);
Ytemp=repmat(y,1,Nterms/numel(y))'; 
SigTemp=z';SigTemp=SigTemp(:);
[xq,yq]=meshgrid(xout,yout);
%A=griddata(Xtemp,Ytemp,SigTemp,xq,yq,'cubic');
F=scatteredInterpolant(Xtemp,Ytemp,SigTemp,'natural','linear');
A = F(xq,yq);
% figure('Visible','on')
% plot(Xtemp,Ytemp,'.');
% hold on
% contourf(xq,yq,A)
% hold off


result=struct('Z',A,'X',xout,'Y',yout);
end
