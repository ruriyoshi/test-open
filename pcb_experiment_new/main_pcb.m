clear all

n = 50;
trange = 400:600;

for i = 1337 %[1333,1336:1340,1342:1345,1346:1348,1350,1352]%[938,940,941,952,953,955]
    load(strcat('C:\Users\uswk0\OneDrive - g.ecc.u-tokyo.ac.jp\data\before_picture\a039_',num2str(i),'.mat'));
   
    %磁気面のみ
    %plot_psi_multi_pcb(data2D,grid2D,false,true,i);
   
    %fitrate_200ch:x点もろもろ&磁気軸のデータ計算コード
    %plot_psi_multi_pcb_for_axis:磁気面 & x点 & 磁気軸描画コード
    [psi_pr, fitrate, xJt, xEt, xeta, xpos,max_psi_left_pos, max_psi_right_pos] = fitrate_200ch(data2D,grid2D,shot,n,trange);
    plot_psi_multi_pcb_for_axis(data2D,grid2D,false,true,i,n,trange,xpos,max_psi_left_pos, max_psi_right_pos);
    saveas(gcf,strcat('C:\Users\uswk0\OneDrive - g.ecc.u-tokyo.ac.jp\data\picture\a039_',num2str(i),'.png'))
   
   
    %plot_merging_rate:合体率抵抗率時間発展描画コード
    plot_merging_rate(psi_pr, fitrate, xJt, xEt, xeta,trange);
    saveas(gcf,strcat('C:\Users\uswk0\OneDrive - g.ecc.u-tokyo.ac.jp\data\picture\a039_',num2str(i),'_fitrate.png'))%
   %}
%%%%


%     helicity_calculation(data2D,rot90(grid2D.rq(:,1),3),grid2D.zq(1,:),rot90(grid2D.rq(:,1),3),grid2D.zq(1,:));
%     lambda_calculation(data2D,grid2D);
%     rogowski(230412,10,300,i);
%     saveas(gcf,strcat('/Users/yunhancai/Downloads/rogo_230225_',num2str(i),'.png'))%,
%    close
%     disp(i)
    clearvars -except i n trange
end

%plot_B_z_in_time(B_z,ch_dist,350,600)
%plot_B_z_vs_z(B_z_new,440,r_bz_probe,z_bz_probe,ch_dist);
%[psi_1,psi_2,psi_common] = merging_rate2(data2D.psi,grid2D.zq(1,:),rot90(grid2D.rq(:,1),3),450-300,600-300);
%plot_psi_at_t(Bzcoil,r_bz_probe,z_bz_probe,1,false,true,true,false)
%plot_B_z_in_rz(B_z_new,r_bz_probe,z_bz_probe);
%plot_B_z_in_rz(B_t_new,r_bt_probe,z_bt_probe);