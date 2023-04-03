clearvars
start=460;
frame=50;
elapsedtime=1:frame;
trange=[1:frame]+start;

IDXlist=[2911:2913 2925 2926 2927 2931 2933];
tstartlist=[472 472 472 477 477 477 482 482];
tendlist=[482 482 482 487 487 487 500 500];
tstartlist=tstartlist-start;
tendlist=tendlist-start;

time=1:20;
fitrate=readmatrix('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\過去データ解析\221229\data.xlsx','Sheet','fitrate_t0','Range','A1:H20');
afit=fitrate(:,1:3);
bfit=fitrate(:,4:6);
cfit=fitrate(:,7:8);
wide=readmatrix('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\過去データ解析\221229\data.xlsx','Sheet','wide_t0','Range','B2:I21');
awide=wide(:,1:3);
bwide=wide(:,4:6);
cwide=wide(:,7:8);
slope=readmatrix('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\過去データ解析\221229\data.xlsx','Sheet','diff_t0','Range','B2:I21');
slope=slope.*180./pi;
aslope=slope(:,1:3);
bslope=slope(:,4:6);
cslope=slope(:,7:8);
eta=readmatrix('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\過去データ解析\221229\data.xlsx','Sheet','eta_t0','Range','B2:I21');
eta=-1.*eta;
aeta=eta(:,1:3);
beta=eta(:,4:6);
ceta=eta(:,7:8);
aeta(17,2)=NaN;


%fitrate
% load('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\220801中間審査\data\fitrate.mat','fitrate','afit','arange','merge_ind');
%fitrate-plot
figure
plot(time,mean(afit(:,1:3),2),'k-+','LineWidth',1)
hold on
plot(time,mean(bfit(:,1:3),2),'r-+','LineWidth',1)
plot(time,mean(cfit(:,1:2),2),'b-+','LineWidth',1)
 ar1 = area(time,[max(afit(:,1:3),[],2), min(afit(:,1:3),[],2)-max(afit(:,1:3),[],2)]);
set(ar1(1),'FaceColor','None','LineStyle','None')
set(ar1(2),'FaceColor','k','FaceAlpha',0.1,'LineStyle','None')
 ar2 = area(time,[max(bfit(:,1:3),[],2), min(bfit(:,1:3),[],2)-max(bfit(:,1:3),[],2)]);
set(ar2(1),'FaceColor','None','LineStyle','None')
set(ar2(2),'FaceColor','r','FaceAlpha',0.1,'LineStyle','None')
 ar3 = area(time,[max(cfit(:,1:2),[],2), min(cfit(:,1:2),[],2)-max(cfit(:,1:2),[],2)]);
set(ar3(1),'FaceColor','None','LineStyle','None')
set(ar3(2),'FaceColor','b','FaceAlpha',0.1,'LineStyle','None')
hold off
ylabel('Merge ratio')
xlabel('Time [us]')
xlim([0 17])%xlim([0 30])
ylim([0 1])
legend('(a)','(b)','(c)','Location','northoutside')
ha1 = gca;
ha1.LineWidth = 1;
ha1.FontSize=13;

%angle-plot
figure
plot(time,mean(aslope(:,1:3),2),'k-+','LineWidth',1)
hold on
plot(time,mean(bslope(:,1:3),2),'r-+','LineWidth',1)
plot(time,mean(cslope(:,1:2),2),'b-+','LineWidth',1)
 ar1 = area(time,[max(aslope(:,1:3),[],2), min(aslope(:,1:3),[],2)-max(aslope(:,1:3),[],2)]);
set(ar1(1),'FaceColor','None','LineStyle','None')
set(ar1(2),'FaceColor','k','FaceAlpha',0.1,'LineStyle','None')
 ar2 = area(time,[max(bslope(:,1:3),[],2), min(bslope(:,1:3),[],2)-max(bslope(:,1:3),[],2)]);
set(ar2(1),'FaceColor','None','LineStyle','None')
set(ar2(2),'FaceColor','r','FaceAlpha',0.1,'LineStyle','None')
 ar3 = area(time,[max(cslope(:,1:2),[],2), min(cslope(:,1:2),[],2)-max(cslope(:,1:2),[],2)]);
set(ar3(1),'FaceColor','None','LineStyle','None')
set(ar3(2),'FaceColor','b','FaceAlpha',0.1,'LineStyle','None')
hold off
ylabel('sheet angle [°]')
xlabel('Time [us]')
xlim([0 17])%xlim([0 30])
ylim([-0.55 0.1].*180./pi)
legend('(a)','(b)','(c)','Location','northoutside')
ha2 = gca;
ha2.LineWidth = 1;
ha2.FontSize=13;

%wide-plot
figure
awide=awide./2;
bwide=bwide./2;
cwide=cwide./2;

plot(time,mean(awide(:,1:3),2),'k-+','LineWidth',1)
hold on
plot(time,mean(bwide(:,1:3),2),'r-+','LineWidth',1)
plot(time,mean(cwide(:,1:2),2),'b-+','LineWidth',1)
 ar1 = area(time,[max(awide(:,1:3),[],2), min(awide(:,1:3),[],2)-max(awide(:,1:3),[],2)]);
set(ar1(1),'FaceColor','None','LineStyle','None')
set(ar1(2),'FaceColor','k','FaceAlpha',0.1,'LineStyle','None')
 ar2 = area(time,[max(bwide(:,1:3),[],2), min(bwide(:,1:3),[],2)-max(bwide(:,1:3),[],2)]);
set(ar2(1),'FaceColor','None','LineStyle','None')
set(ar2(2),'FaceColor','r','FaceAlpha',0.1,'LineStyle','None')
 ar3 = area(time,[max(cwide(:,1:2),[],2), min(cwide(:,1:2),[],2)-max(cwide(:,1:2),[],2)]);
set(ar3(1),'FaceColor','None','LineStyle','None')
set(ar3(2),'FaceColor','b','FaceAlpha',0.1,'LineStyle','None')
hold off
ylabel('delta [m]')
xlabel('Time [us]')
xlim([0 17])%xlim([0 30])
ylim([0.02 0.032])%ylim([0.045 0.065])%ylim([0.01 0.045])
legend('(a)','(b)','(c)','Location','northoutside')
ha1 = gca;
ha1.LineWidth = 1;
ha1.FontSize=13;
%eta-plot
figure
plot(time,mean(aeta(:,1:3),2).*1e+3,'k-+','LineWidth',1)
hold on
plot(time,mean(beta(:,1:3),2).*1e+3,'r-+','LineWidth',1)
plot(time,mean(ceta(:,1:2),2).*1e+3,'b-+','LineWidth',1)
 ar1 = area(time,[max(aeta(:,1:3),[],2).*1e+3, min(aeta(:,1:3),[],2).*1e+3-max(aeta(:,1:3),[],2).*1e+3]);
set(ar1(1),'FaceColor','None','LineStyle','None')
set(ar1(2),'FaceColor','k','FaceAlpha',0.1,'LineStyle','None')
 ar2 = area(time,[max(beta(:,1:3),[],2).*1e+3, min(beta(:,1:3),[],2).*1e+3-max(beta(:,1:3),[],2).*1e+3]);
set(ar2(1),'FaceColor','None','LineStyle','None')
set(ar2(2),'FaceColor','r','FaceAlpha',0.1,'LineStyle','None')
 ar3 = area(time,[max(ceta(:,1:2),[],2).*1e+3, min(ceta(:,1:2),[],2).*1e+3-max(ceta(:,1:2),[],2).*1e+3]);
set(ar3(1),'FaceColor','None','LineStyle','None')
set(ar3(2),'FaceColor','b','FaceAlpha',0.1,'LineStyle','None')
hold off
ylabel('resistivity [mΩ・m]')
xlabel('Time [us]')
xlim([0 17])%xlim([0 30])
ylim([0 2])
legend('(a)','(b)','(c)','Location','northoutside')
ha1 = gca;
ha1.LineWidth = 1;
ha1.FontSize=13;

%IDX2911
data_a1=readmatrix('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\過去データ解析\221229\data.xlsx','Sheet','plot2911');
time=data_a1(:,2);
fit_a1=data_a1(:,3);
Jt_a1=-1.*data_a1(:,4).*1e-6;%[MA/m^{2}]
Et_a1=-1.*data_a1(:,5).*1e-3;%[kV/m]
eta_a1=data_a1(:,6).*1e+3;%[mΩ・m]
slope_a1=data_a1(:,7).*180./pi;%[degree]
wide_a1=data_a1(:,8)./2;

%(a)2911のみ
figure
% f1=figure;
% f1.WindowState='maximized';
subplot(5,1,1)
plot(time,fit_a1,'k-+','LineWidth',1)
ylabel('Merge ratio')
xlim([0 17])%xlim([0 30])
ylim([0 1])
xticklabels({})
ha1 = gca;
ha1.LineWidth = 1;
ha1.FontSize=13;

Jt_a1([1 2])=NaN;
Et_a1([1 2])=NaN;
eta_a1([1 2])=NaN;
Jt_a1(13:end)=NaN;
Et_a1(13:end)=NaN;
eta_a1(13:end)=NaN;
subplot(5,1,2)
yyaxis left
plot(time,Jt_a1,'b-+','LineWidth',1)
ylabel('Jt [MA/m^{2}]')
xlim([0 20])
ylim([0 1.1])%ylim([0 1.5])%pcb
xticklabels({})
yyaxis right
plot(time,Et_a1,'r-+','LineWidth',1)
ylabel('Et [kV/m]')
xlim([0 17])
ylim([0 0.35])%ylim([0 0.4])%pcb
xticklabels({})
ha2 = gca;
ha2.LineWidth = 1;
ha2.FontSize=13;

subplot(5,1,3)
plot(time,eta_a1,'k-+','LineWidth',1)
% hold on
% yline(0.1,'b--')%古典抵抗
% hold off
ylabel('η [mΩ・m]')
% xlabel('Time [us]')
xlim([0 17])%xlim([0 30])
ylim([0 1.2])%ylim([0 0.6])%pcb
xticklabels({})
ha3 = gca;
ha3.LineWidth = 1;
ha3.FontSize=13;

slope_a1(15:end)=NaN;
subplot(5,1,4)
plot(time,slope_a1,'k-+','LineWidth',1)
ylabel('angle [°]')
xlim([0 17])%xlim([0 30])
ylim([-25 0])%([-0.55 0].*180./pi)
xticklabels({})
ha4 = gca;
ha4.LineWidth = 1;
ha4.FontSize=13;

wide_a1(13:end)=NaN;
subplot(5,1,5)
plot(time,wide_a1,'k-+','LineWidth',1)
ylabel('δ [m]')
xlabel('Time [us]')
xlim([0 17])%xlim([0 30])
ylim([0.02 0.026])%ylim([0.01 0.045])
ha5 = gca;
ha5.LineWidth = 1;
ha5.FontSize=13;