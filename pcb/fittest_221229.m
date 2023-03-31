%IDX2911,t=473us, 補間方法を変化させた場合のBz生信号
load('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\過去データ解析\221229\vq_t473_fit_lowess.mat');
vq_fit=vq;
clear vq
load('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\過去データ解析\221229\vq_t473_rbfinterp.mat');
vq_rbf=vq;
clear vq

load('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\過去データ解析\221229\data2Dcalc_fit_lowess.mat','grid2D');
load('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\過去データ解析\221229\data2Dcalc_fit_lowess.mat','trange');
load('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\過去データ解析\221229\data2Dcalc_fit_lowess.mat','data2D');
data2D_fit=data2D;
clear data2D
load('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\過去データ解析\221229\data2Dcalc_rbfinterp.mat','data2D');
data2D_rbf=data2D;
clear data2D

% figure
% contourf(grid2D.zq(1,:),grid2D.rq(:,1),vq_fit,30,'LineStyle','none')
% colormap(jet)
% axis image
% axis tight manual
% caxis([-0.1,0.1])
% colorbar('Location','eastoutside')
% hold on
% contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D_fit.psi(:,:,4)),[-20e-3:0.2e-3:40e-3],'black')
% hold off 
% title('vq at t=473us by fit lowess')
% 
% figure
% contourf(grid2D.zq(1,:),grid2D.rq(:,1),vq_rbf,30,'LineStyle','none')
% colormap(jet)
% axis image
% axis tight manual
% caxis([-0.1,0.1])
% colorbar('Location','eastoutside')
% hold on
% contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D_rbf.psi(:,:,4)),[-20e-3:0.2e-3:40e-3],'black')
% hold off 
% title('vq at t=473us by RBF')

figure
contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D_fit.Jt(:,:,4),30,'LineStyle','none')
colormap(jet)
axis image
axis tight manual
caxis([-3e+6,3e+6])
colorbar('Location','eastoutside')
hold on
contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D_fit.psi(:,:,4)),[-20e-3:0.2e-3:40e-3],'black')
hold off 
title('Jt at t=473us by fit lowess')

figure
contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D_rbf.Jt(:,:,4),30,'LineStyle','none')
colormap(jet)
axis image
axis tight manual
caxis([-3e+6,3e+6])
colorbar('Location','eastoutside')
hold on
contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D_fit.psi(:,:,4)),[-20e-3:0.2e-3:40e-3],'black')
hold off 
title('Jt at t=473us by RBF')

