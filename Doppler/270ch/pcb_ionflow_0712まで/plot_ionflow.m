function [] = plot_ionflow(V_i,absV,T_i,date,expval,ICCD,pathname,factor,mpoints,overlay_plot,save_fig,cal_type)
%V_i/absV/T_i/������/��������/ICCD�ϐ�/pathname/���T�C�Y(���l:0.05�Ȃ�)/�v���_�z��/���C�ʂɏd�˂�/fig��ۑ�/�����v�Z���@

time = round(ICCD.trg+ICCD.exp_w/2);%�v������
if overlay_plot%���C�ʂƏd�˂�
    unit = 1E-2;%�P�ʂ�cm����m�ɕϊ�
else
    unit = 1;%�P�ʕϊ����Ȃ�
    figure('Position',[600 150 300 600])
end
mpoints.z = mpoints.z * unit;
mpoints.r = mpoints.r * unit;
factor = factor * unit;

% %���x�J���[�}�b�v
% T_icon = repmat(T_i(:,1),1,5);%�������}�p���σC�I�����x
% zcon = mpoints.z(1,1);
% s = pcolor([zcon-2*unit zcon-1*unit zcon zcon+1*unit zcon+2*unit],mpoints.r,T_icon);
% s.FaceColor = 'interp';
% s.EdgeAlpha = 0;
% colormap('jet')
% colorbar
% clim([0 150])
% hold on

r_start = 1;
r_end = 5;
% plot(mpoints.z,mpoints.r,'xr');%�h�b�v���[�v���[�u�v���_��\��
plot(mpoints.z(r_start:r_end,1),mpoints.r(r_start:r_end,1),'xk');%�h�b�v���[�v���[�u�v���_��\��
hold on
if mpoints.n_z == 1
    q = quiver(mpoints.z(r_start:r_end,1),mpoints.r(r_start:r_end,1),V_i(r_start:r_end,1)*factor,V_i(r_start:r_end,2)*factor);
end
if mpoints.n_z == 2
    q = quiver(mpoints.z(r_start:r_end,1),mpoints.r(r_start:r_end,1),[V_i(r_start:r_end,1),V_i(r_start:r_end,3)]*factor,[V_i(r_start:r_end,2),V_i(r_start:r_end,4)]*factor);
end
q.LineWidth = 2;
q.MaxHeadSize = 10;
q.AutoScale = 'off';
q.Color = 'r';
absV = round(absV,1);
for j = 1:mpoints.n_z
    for i = r_start:r_end
        txt = text(mpoints.z(i,j)+0.8*unit,mpoints.r(i,j)-0.2*unit,num2str(absV(i,j)));
        txt.FontSize = 18;
        txt.Color = 'r';
        txt.FontWeight = 'bold';
    end
end
if overlay_plot
    s.FaceAlpha = 0.3;%�����x(0�Ŋ��S�ɓ���)
    xlim([-0.17 0.17])
    ylim([0.06 0.33])
    % xlim([-0.04 0.02])
    % ylim([0.125 0.275])
    % xticks([-0.04 -0.01 0.02])
    % yticks(0.15:0.05:0.25)
    daspect([1 1 1])
    ax = gca;
    ax.FontSize = 16;
    title([num2str(time),' [us]'])
    % ax.FontSize = 18;
    % title(['shot', num2str(ICCD.shot),'-',num2str(time),' [us] Ion Flow [km/s]',newline,'PF1 = ',num2str(expval.PF1), ...
    %     ' [kV], PF2 = ',num2str(expval.PF2),' [kV], TF = ',num2str(expval.TF),' [kV], EF = ',num2str(expval.EF),' [A]'] ...
    %     ,'Color','black','FontWeight','bold')
    if save_fig
        if not(exist([pathname.fig,'/',cal_type,'_psi/',num2str(date)],'dir'))
            mkdir(sprintf("%s/%s_psi",pathname.fig,cal_type), sprintf("%s", num2str(date)));
        end
        saveas(gcf,[pathname.fig,'/',cal_type,'_psi/',num2str(date),'/','shot', num2str(ICCD.shot),'_',num2str(time),'us_PF1_',num2str(expval.PF1), ...
        'kV_PF2_',num2str(expval.PF2),'kV_TF_',num2str(expval.TF),'kV_EF_',num2str(expval.EF),'A.png'])
        hold off
        close
    else
        hold off
    end
else
    xlim([min(mpoints.z,[],'all')-2.5 max(mpoints.z,[],'all')+2.5])
    ylim([min(mpoints.r,[],'all')-2.5 max(mpoints.r,[],'all')+2.5])
    % title(['shot', num2str(ICCD.shot),'-',num2str(time),' [us]', newline, 'Ion Flow [km/s]'],'Color','black','FontWeight','bold')
    title(['shot', num2str(ICCD.shot),'-',num2str(time),' [us] Ion Flow [km/s]',newline,'PF1 = ',num2str(expval.PF1), ...
    ' [kV], PF2 = ',num2str(expval.PF2),' [kV]',newline,'TF = ',num2str(expval.TF),' [kV], EF = ',num2str(expval.EF),' [A]'] ...
    ,'Color','black','FontWeight','bold')
    xlabel('Z [cm]')
    ylabel('R [cm]')
    ax = gca;
    ax.FontSize = 10;
    grid on
    daspect([1 1 1])
    if save_fig
        if not(exist([pathname.fig,'/',cal_type,'/',num2str(date)],'dir'))
            mkdir(sprintf("%s/%s", pathname.fig,cal_type), sprintf("%s", num2str(date)));
        end
        saveas(gcf,[pathname.fig,'/',cal_type,'/',num2str(date),'/','shot', num2str(ICCD.shot),'_',num2str(time),'us_PF1_',num2str(expval.PF1), ...
        'kV_PF2_',num2str(expval.PF2),'kV_TF_',num2str(expval.TF),'kV_EF_',num2str(expval.EF),'A.png'])
        hold off
        close
    else
        hold off
    end
end
end