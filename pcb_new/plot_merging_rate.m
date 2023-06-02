
function plot_merging_rate(psi_pr, fitrate, xJt, xEt, xeta,trange)

%

%合体率のplot
figure('Position',[700,350,500,300],'visible','on')
yyaxis left
plot(trange,psi_pr(1,:),'k--+','LineWidth',1)
hold on
plot(trange,psi_pr(2,:),'k--+','LineWidth',1)
plot(trange,psi_pr(3,:),'b--+','LineWidth',1)
hold off
ylabel('Psi [Wb]')
xlabel('Time [us]')
ylim([0 inf])
xlim([437 473])
yyaxis right
plot(trange,fitrate,'r-+','LineWidth',1)
ylabel('Fitrate')
legend('Psi_{private1}','Psi_{private2}','Psi_{common}','Fitrate','Location','eastoutside')
ylim([0 1])
xlim([456 475])
ha1 = gca;
ha1.LineWidth = 1;
ha1.FontSize=10;
% saveas(gcf,strcat('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\230123\fitrate\',num2str(shot),'_fitrate_TF',num2str(TF),'kV.png'))
% save(strcat('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\230123\fitrate\fitrate_',num2str(shot),'.mat'),'date','i_EF','fitrate','n','psi_pr','shot','TF','tfshot','trange')
% close


%合体率とX点Jt,Et,etaのplot
figure('Position',[0,0,300,700],'visible','on')
subplot(3,1,1)
plot(trange,fitrate,'k-+','LineWidth',1)
ylabel('Merging Ratio')
xlabel('Time [us]')
ylim([0 1])
xlim([456 472])
ha1 = gca;
ha1.LineWidth = 1;
ha1.FontSize=10;

subplot(3,1,2)
yyaxis left
plot(trange,-1.*xJt,'b-+','LineWidth',1)
ylabel('Jt [A/m^{2}]')
xlabel('Time [us]')
xlim([456 472])
ylim([0 inf])
yyaxis right
plot(trange,-1.*xEt,'r-+','LineWidth',1)
ylabel('Et [V/m]')
xlim([455 472])
ylim([0 inf])
ha2 = gca;
ha2.LineWidth = 1;
ha2.FontSize=10;

subplot(3,1,3)
plot(trange,xeta,'k-+','LineWidth',1)
ylabel('η [Ω m]')
xlabel('Time [us]')
xlim([456 472])
ylim([0 6e-3])
ha1 = gca;
ha1.LineWidth = 1;
ha1.FontSize=10;

% saveas(gcf,strcat('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\230123\xpoint\',num2str(shot),'_xpoint_TF',num2str(TF),'kV.png'))
% save(strcat('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\230123\xpoint\xpoint_',num2str(shot),'.mat'),'date','i_EF','fitrate','n','psi_pr','shot','TF','tfshot','trange','xEt','xeta','xJt','xpos')
% close

end