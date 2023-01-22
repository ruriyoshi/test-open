function doppler=doppler288ch(filename,date)
tmax=300
EM_max=1.e8
lambda0=480.6%656.3%471.3%656.3%480.6%587.562%486.133%656.3%486.133%468.57%486.133%nm
mass=39.95%4.%12.01%39.95%1.
resolution=-1.420387e-12*lambda0^3 - 2.156031e-09*lambda0^2 + 1.250038e-06*lambda0 + 3.830769e-03%0.0037714% -0.000000000001420387*lambda0^3 - 0.000000002156031*lambda0^2 + 0.000001250038*lambda0 + 0.003830769

z=readmatrix("z.txt")*1e-3;
%[file,path] = uigetfile('*.asc')
%filename=fullfile(path,file)

[d,delimiterOut]=importdata(filename);d=d(:,2:end)';
%[file,path] = uigetfile('*.asc')
%filenamebg=fullfile(path,file)
filenamebg=strcat('I:\Doppler\288CH\20',num2str(date),'\bg.asc');
bg=importdata(filenamebg,delimiterOut);bg=bg(:,2:end)';d=d-bg;
if date==211217
calibfile="Ar_calibration_20211217.txt";
elseif date==211219
    calibfile="H_calibration_20211219.txt";
    lambda0=486.133;%656.3%471.3%656.3%480.6%587.562%486.133%656.3%486.133%468.57%486.133%nm
mass=1;%4.%12.01%39.95%1.
resolution=-1.420387e-12*lambda0^3 - 2.156031e-09*lambda0^2 + 1.250038e-06*lambda0 + 3.830769e-03%0.0037714% -0.000000000001420387*lambda0^3 - 0.000000002156031*lambda0^2 + 0.000001250038*lambda0 + 0.003830769
elseif date==211223 ||  date==211224
     calibfile="Ar_calibration_20211223-24.txt";
end     
calib=readmatrix(calibfile);
%calib=readmatrix("Ar_calibration_20211223-24.txt")
CH=squeeze(calib(:,1));
Center=squeeze(calib(:,2));
Smile=squeeze(calib(:,3))+8;
relative=squeeze(calib(:,5));
instru=squeeze(calib(:,6));



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
 lambdaB=lambdaA ;
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
    passive.Ti(i)=1.69e8*mass*(resolution*c.c1/sqrt(2)*sqrt(2.*log(2.))/lambda0)^2-Ti_instru(i);
    passive.Timax(i)=1.69e8*mass*(resolution*conf(2,3)/sqrt(2)*sqrt(2.*log(2.))/lambda0)^2-Ti_instru(i);
    passive.Timin(i)=1.69e8*mass*(resolution*conf(1,3)/sqrt(2)*sqrt(2.*log(2.))/lambda0)^2-Ti_instru(i);
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
%xline(CH(separation),':')
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

num=16;
%yy=[1:num]/num*(edge-min(p))+min(p);
yy=[0:num]/(num)*(edge-0.06)+0.06;
spectra_interp=zeros([numel(yy),numel(lambdaA),numel(z)]);
%grid2D=meshgrid(lambdaA,yy);
for i=1:numel(z) 
    [Yin,I]= sort(pall(CH(zall(CH)==zall(1,i))));
    Xin=lambda(:,zall(CH)==zall(1,i));
    Zin=spectra(:,zall(CH)==zall(1,i));
    Yin=[Yin', edge];
    Xin=[Xin(:,I),lambdaB'];
    Zin=[Zin(:,I),zeros(numel(lambdaB),1)];
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

 
 % spectra_interp=smoothdata(spectra_interp,1,'movmean',num/16);
% spectra_interp=smoothdata(spectra_interp,2,'movmean',5);
% spectra_interp=smoothdata(spectra_interp,3,'movmean',1);
% dy=yy(2)-yy(1);
reconst_int=spectra_interp*0;
[~,dndy,~]=gradient(spectra_interp,1,yy,1);
dy=diff(yy);
for i=1:numel(z)
 for j=1:num+1
  for k=j+1:num+1 
     reconst_int(j,:,i)=reconst_int(j,:,i)-1./pi.*dndy(k,:,i)./sqrt(yy(k).^2-yy(j).^2).*dy(k-1);
  end
 end
end
% reconst_int=smooth(reconst_int,[5,num/16,1])
% 
 emission=zeros([numel(yy),numel(z)]);
 Ti.Ti2D=zeros([numel(yy),numel(z)]);
 Ti.max=zeros([numel(yy),numel(z)]);
 Ti.min=zeros([numel(yy),numel(z)]);
 Ti_instru2=Ti_instru*0;
 for i=1:numel(z)
     Ti_instru2(zall(CH)==zall(1,i))=sum(Ti_instru(zall(CH)==zall(1,i)))/sum(zall(CH)==zall(1,i));
 end 
 for j=1:numel(z)
% figure('Visible','on')
    for i=1:numel(yy)  

%    if mod(i, 36)==1 
%        figure('Visible','on')
%    end
 % subplot(ceil(sqrt(numel(yy))),ceil(numel(yy)/sqrt(numel(yy))),i)
    input=squeeze(reconst_int(i,:,j));%input(input<=0)=-input(input<=0);
    c=fit(lambdaA',input','gauss1','StartPoint',[max(input),480.6,0.002]);   
    conf=confint(c);
    Ti.Ti2D(i,j)=1.69e8*mass*(2*c.c1/sqrt(2)*sqrt(2.*log(2.))/lambda0)^2-Ti_instru2(j);
    Ti.max(i,j)=1.69e8*mass*(2*conf(2,3)/sqrt(2)*sqrt(2.*log(2.))/lambda0)^2-Ti_instru2(j);
    Ti.min(i,j)=1.69e8*mass*(2*conf(1,3)/sqrt(2)*sqrt(2.*log(2.))/lambda0)^2-Ti_instru2(j);
    emission(i,j)=resolution*sum(input);
    checker=(abs(c.b1-lambda0) < 0.3) & (c.a1  >= 1e5) %&(abs(Ti.max(i,j)-Ti.min(i,j)) < Ti.Ti2D(i,j))
    %checker=(abs(c.b1-lambda0) < 0.3) & & %*float(abs(Ti.max(i,j)-Ti.min(i,j)) lt Ti.Ti2D(i,j)+Ti_instru2[j])%*(emission_local(i,j) gt EM_max*0.1);
 if checker==false
     Ti.Ti2D(i,j)=NaN;
     Ti.max(i,j)=NaN;
     Ti.min(i,j)=NaN;
     emission(i,j)=NaN;
%    plot(c,lambdaA',input',':')
 else 
%    plot(c,lambdaA',input')
 end
 end 
 end  
 doppler=struct('ti_2d',Ti.Ti2D, 'emission',emission, 'z',z,'yy',yy);
 close all
end