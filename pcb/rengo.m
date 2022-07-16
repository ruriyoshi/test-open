% maxslope=[-0.5036 -0.5490 -0.5160 -0.1692 -0.1775 -0.1321 +0.0619 +0.0826];
% maxslope=atan(maxslope).*180./pi;
% 
% figure
% plot(1,mean(maxslope(1:3)),'*')
% hold on
% plot(1.1,mean(maxslope(4:6)),'*')
% plot(1.2,mean(maxslope(7:8)),'*')
% yline(0,'k--')
% hold off
% legend('No Cr','Cr 454 us','Cr 450 us')
% xlabel('')
% xlim([0 3])
% ylabel('current sheet slope [degrees]')
% ha1 = gca;
% ha1.LineWidth = 1;
% ha1.FontSize=13;

load('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\過去データ解析\rengo1.mat')
%rengo_xpoint_plot.m-IDX2911-465:510usのデータ
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
plot(time,-1.*xJt.*1e-6,'-o','Color','b')
%yline(0,'k-')
xline(472,'b-',"LineWidth",1)
ylabel('Jt at Xpoint [MA/m^{2}]','FontSize',13,'Color','b','LineWidth',1)
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
