close all

%�ePC�̃p�X���`
run define_path.m

%------�yinput�z-------
date = 230524;%�yinput�z������230524,230526
shotlist_cal = 15:24;%�yinput�z���C��&�t���[�v�Z�V���b�g���X�g15:24,[2:5,7:12]
int_r = 2.5;%�yinput�z�h�b�v���[�v���[�u�v���_r�����Ԋu[cm](2.5)
int_z = 4.2;%�yinput�z�h�b�v���[�v���[�u�v���_z�����Ԋu[cm](4.2)
IDSP.line = 'Ar';%�yinput�z�h�b�v���[�������C��('Ar')
n_CH = 20;%�yinput�z�h�b�v���[�v���[�u�t�@�C�o�[CH��(28)
n_z = 1;%�yinput�z�h�b�v���[�v���[�uz�����f�[�^��(���l)(1)
start_r = 1;%�yinput�z�v���b�g�n��CH
end_r = 5;%�yinput�z�v���b�g�I���CH
trange = 430:590;%�yinput�z���C�v���[�u�v�Z���Ԕ͈�(430:590)
dtacq.num = 39;

mu0 = 4*pi*1e-7;%�^��̓�����
[~,n_shot] = size(shotlist_cal);%�v�Z�V���b�g��

%�������O�ǂݎ��
[exp_log,index,begin_row,end_row] = load_log(date);
if isempty(begin_row)
    return
end

%�S�V���b�g�̃f�[�^������z�������
V_i_all = zeros(n_CH/4,2*n_z,n_shot);
absV_all = zeros(n_CH/4,n_z,n_shot);
T_i_all = zeros(n_CH/4,n_z,n_shot);
magpres_t_all = zeros(100,1,n_shot);
magpres_tp_all = zeros(100,1,n_shot);
magpres_z_all = zeros(100,1,n_shot);
magpres_all = zeros(100,1,n_shot);
diff_magpres_all = zeros(100-1,1,n_shot);

for i = 1:n_shot
    i_log = shotlist_cal(1,i) + begin_row - 1;
    IDSP.shot = exp_log(i_log,index.shot);%�V���b�g�ԍ�
    a039shot = exp_log(i_log,index.a039);%a039�V���b�g�ԍ�
    a039tfshot = exp_log(i_log,index.a039_TF);%a039TF�V���b�g�ԍ�
    IDSP.trg = exp_log(i_log,index.IDSP_trg);%IDSP�g���K����
    IDSP.exp_w = exp_log(i_log,index.IDSP_exp_w);%IDSP�I������
    IDSP.gain = exp_log(i_log,index.IDSP_gain);%Andor gain
    time = round(IDSP.trg+IDSP.exp_w/2);%�v������
    min_r = exp_log(i_log,index.minR);%�h�b�v���[�v���[�u�v���_�ŏ�r���W
    min_z = exp_log(i_log,index.minZ);%�h�b�v���[�v���[�u�v���_�ŏ�z���W
    cut_z_magpres = min_z;%���C���v�Zz���W
    if dtacq.num == 39
        dtacq.shot = a039shot;
        dtacq.tfshot = a039tfshot;
    end
    %�h�b�v���[�v���[�u�v���_�z��𐶐�
    mpoints = make_mpoints(n_CH,min_r,int_r,n_z,min_z,int_z);
    %�ۑ��ς݃C�I�����x�A�t���[��ǂݎ��
    % [V_i_all(:,:,i),absV_all(:,:,i),T_i_all(:,:,i)] = load_ionflow(date,IDSP,pathname);
    [V_i_all(:,:,i),absV_all(:,:,i),T_i_all(:,:,i),~,~,~,~,~,~,~,~] = load_ionvdist(date,IDSP,pathname);
    %�ۑ��ςݎ���f�[�^��ǂݎ��
    [grid2D,data2D,ok_z,ok_r] = load_pcb200ch(date,dtacq,pathname);
    [~,data2D_tfoffset,~,~] = load_pcb200ch_with_tfoffset(date,dtacq,pathname);
    %���C���v�Z
    cut_z_magpres = cut_z_magpres*1e-2;%�P�ʂ�[cm]����[m]�ɕϊ�
    i_pcb = time-trange(1)+1; %����f�[�^�̎���index
    IDX = knnsearch(grid2D.zq(1,:).',cut_z_magpres);%grid2D.zq�̒��ōł�cut_z�ɋ߂��Z���ԍ����擾
    z = grid2D.zq(1,IDX)*1e2;%z���W[cm]
    [nR,~,~] = size(data2D.Bz);
    magpres_z = zeros(nR,1);
    magpres_t = zeros(nR,1);
    magpres_tp = zeros(nR,1);
    magpres = zeros(nR,1);
    for j = 1:nR
        magpres_z(j,1) = data2D.Bz(j,IDX,i_pcb)^2/(2*mu0);
        magpres_t(j,1) = data2D_tfoffset.Bt(j,IDX,i_pcb)^2/(2*mu0);
        magpres_tp(j,1) = data2D.Bt(j,IDX,i_pcb)^2/(2*mu0);
        magpres(j,1) = magpres_z(j,1) + magpres_t(j,1);
    end
    magpres_t_all(:,1,i) = magpres_t;
    magpres_tp_all(:,1,i) = magpres_tp;
    magpres_z_all(:,1,i) = magpres_z;
    magpres_all(:,1,i) = magpres;

    for k = 1:99
        dr = grid2D.rq(2,1) - grid2D.rq(1,1);
        diff_magpres_all(k,1,i) = (magpres_all(k+1,1,i) - magpres_all(k,1,i))/(grid2D.rq(k+1,1) - grid2D.rq(k,1));
    end
end

if not(exist([pathname.mat,'/ionflow_magpres/',num2str(date)],'dir'))
    mkdir(sprintf("%s/ionflow_magpres", pathname.mat), sprintf("%s", num2str(date)));
end
save([pathname.mat,'/ionflow_magpres/',num2str(date),'/shot',num2str(shotlist_cal(1)),'-',num2str(shotlist_cal(n_shot)),'.mat'],'grid2D','mpoints','V_i_all','T_i_all','magpres_t_all','magpres_tp_all','magpres_z_all','magpres_all','diff_magpres_all')

% V_i_mean = mean(V_i_all,3);
% V_i_sigma = std(V_i_all,0,3);
% magpres_t_mean = mean(magpres_t_all,3);
% magpres_t_sigma = std(magpres_t_all,0,3);
% magpres_z_mean = mean(magpres_z_all,3);
% magpres_z_sigma = std(magpres_z_all,0,3);
% magpres_mean = mean(magpres_all,3);
% magpres_sigma = std(magpres_all,0,3);
% figure('Position',[600 800 600 600])
% %���c��
% yyaxis left
% errorbar(mpoints.r(start_r:end_r),-V_i_mean(start_r:end_r,2),V_i_sigma(start_r:end_r,2),'LineWidth',2)%Vr
% xlabel('R [cm]')
% ylabel('Ion Velocity [km/s]')
% left_y_upper = max(-V_i_mean(start_r:end_r,2) + V_i_sigma(start_r:end_r,2));
% left_y_lower = min(-V_i_mean(start_r:end_r,2) - V_i_sigma(start_r:end_r,2));
% left_ylim = round(max([left_y_upper -left_y_lower]))+1;
% ylim([-left_ylim left_ylim])
% yline(0,'--','LineWidth',3);
% %�E�c��
% yyaxis right
% unit = 1E2;
% errorbar(grid2D.rq(:,1)*unit,magpres_mean(:,1),magpres_sigma(:,1),'LineWidth',2)%���C��
% ylabel('Magneic Pressure [Pa]')
% right_y_upper = max(magpres_mean(:,1) + magpres_sigma(:,1));
% right_y_lower = min(magpres_mean(:,1) - magpres_sigma(:,1));
% right_ylim = round(max([right_y_upper -right_y_lower]))+1;
% ylim([-right_ylim right_ylim])
% xlim([8 20])
% % title([num2str(time),'us'],'Color','black','FontWeight','bold')
% % legend('V_z','V_r')
% ax = gca;
% ax.FontSize = 16;
