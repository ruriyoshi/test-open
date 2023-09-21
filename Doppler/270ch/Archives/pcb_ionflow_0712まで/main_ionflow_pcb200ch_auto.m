%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�V���b�g�ԍ��A�B�e�p�����[�^�Ȃǂ��������O���玩���擾����
%�h�b�v���[�v���[�u�ɂ��C�I�����x�A�t���[�Ƃ��̏u�Ԃ̎��C�ʂ��v���b�g
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%�ePC�̃p�X���`
run define_path.m

%------�yinput�z-------
date = 230524;%�yinput�z������
begin_cal = 15;%�yinput�z���C��&�t���[�v�Z�n��shot�ԍ�(�������OD��)
end_cal = 24;%�yinput�z���C��&�t���[�v�Z�I���shot�ԍ�(�������OD��)(0�ɂ����begin_cal�ȍ~�̓����̑Sshot�v�Z)
int_r = 2.5;%�yinput�z�h�b�v���[�v���[�u�v���_r�����Ԋu[cm](2.5)
int_z = 4.2;%�yinput�z�h�b�v���[�v���[�u�v���_z�����Ԋu[cm](4.2)
ICCD.line = 'Ar';%�yinput�z�h�b�v���[�������C��('Ar')
n_CH = 28;%�yinput�z�h�b�v���[�v���[�u�t�@�C�o�[CH��(28)
n_z = 1;%�yinput�z�h�b�v���[�v���[�uz�����f�[�^��(���l)(1)
%------�ڍאݒ�yinput�z-------
cal_flow = true;%�yinput�z�������v�Z(true,false)
save_flow = true;%�yinput�z�����f�[�^��ۑ�(true,false)
load_flow = false;%�yinput�z�����f�[�^��ǂݍ���(true,false)

cal_pcb = false;%�yinput�z������v�Z(true,false)
save_pcb = false;%�yinput�z����f�[�^��ۑ�(true,false)
load_pcb = true;%�yinput�z����f�[�^��ǂݍ���(true,false)

plot_fit = true;%�yinput�z�K�E�X�t�B�b�e�B���O��\��(true,false)
plot_flow = true;%�yinput�z�������v���b�g(true,false)
plot_psi = true;%�yinput�z���C�ʂ��v���b�g(true,false)
overlay_plot = true;%�yinput�z�����Ǝ��C�ʂ��d�˂�(true,false)

save_fit = true;%�yinput�z�K�E�X�t�B�b�e�B���Opng��ۑ�(true,false)
save_fig = true;%�yinput�z����png��ۑ�(true,false)

show_offset = true;%�yinput�z����offset��\��(true,false)
factor = 0.1;%�yinput�z�C�I���t���[���T�C�Y(���l:0.1�Ȃ�)
dtacq.num = 39;%�yinput�z���C�v���[�udtacq�ԍ�(39)
mesh_rz = 100;%�yinput�z���C�v���[�urz�����̃��b�V����(50)
trange = 430:590;%�yinput�z���C�v���[�u�v�Z���Ԕ͈�(430:590)

%�������O�ǂݎ��
[exp_log,index,begin_row,end_row] = load_log(date);
if isempty(begin_row)
    return
end

%--------���C��&�t���[���v�Z------
start_i = begin_row + begin_cal - 1;
if start_i <= end_row
    if end_cal == 0
        end_i = end_row;%begin_cal�ȍ~�S���v�Z
    elseif end_cal < begin_cal
        error('end_cal must <= begin_cal.')
    elseif begin_row + end_cal - 1 <= end_row
        end_i = begin_row + end_cal - 1;%begin_cal����end_cal�܂Ōv�Z
    else
        error('end_cal must <= %d.', exp_log(end_row,4))
    end
    for i = start_i:end_i
        IDSP.shot = exp_log(i,4);%�V���b�g�ԍ�
        a039shot = exp_log(i,index.a039);%a039�V���b�g�ԍ�
        a039tfshot = exp_log(i,index.a039_TF);%a039TF�V���b�g�ԍ�
        expval.PF1 = exp_log(i,index.PF1);%PF1�d��(kv)
        expval.PF2 = exp_log(i,index.PF2);%PF2�d��(kv)
        expval.TF = exp_log(i,index.TF);%PF2�d��(kv)
        expval.EF = exp_log(i,index.EF);%EF�d��
        IDSP.trg = exp_log(i,index.IDSP_trg);%IDSP�g���K����
        IDSP.exp_w = exp_log(i,index.IDSP_exp_w);%IDSP�I������
        IDSP.gain = exp_log(i,index.IDSP_gain);%Andor gain
        time = round(IDSP.trg+IDSP.exp_w/2);%���C�ʃv���b�g����
        min_r = exp_log(i,index.IDSP_minR);%�yinput�z�h�b�v���[�v���[�u�v���_�ŏ�r���W
        min_z = exp_log(i,index.IDSP_minZ);%�yinput�z�h�b�v���[�v���[�u�v���_�ŏ�z���W
        %�h�b�v���[�v���[�u�v���_�z��𐶐�
        mpoints = make_mpoints(n_CH,min_r,int_r,n_z,min_z,int_z);

        if dtacq.num == 39
            dtacq.shot = a039shot;
            dtacq.tfshot = a039tfshot;
        end
        if cal_flow
            %�C�I�����x�A�t���[���v�Z
            [V_i,absV,T_i] = cal_ionflow(date,ICCD,mpoints,pathname,show_offset,plot_fit,save_fit,save_flow);
        elseif load_flow
            %�ۑ��ς݃C�I�����x�A�t���[��ǂݎ��
            [V_i,absV,T_i] = load_ionflow(date,ICCD,pathname);
        end
        %������v�Z
        if cal_pcb
            [grid2D,data2D,ok_z,ok_r] = cal_pcb200ch(date,dtacq,pathname,mesh_rz,expval,trange,save_pcb);
        elseif load_pcb
            [grid2D,data2D,ok_z,ok_r] = load_pcb200ch(date,dtacq,pathname);
        end
        if not(isempty(data2D))
            %���C�ʂ��v���b�g
            if plot_psi
                plot_psi200ch_at_t(time,trange,grid2D,data2D,ok_z,ok_r);
            end
        end
        %�C�I�����x�A�t���[���v���b�g
        if plot_flow
            if not(isempty(V_i))
                if plot_psi
                    plot_ionflow(V_i,absV,T_i,date,expval,ICCD,pathname,factor,mpoints,overlay_plot,save_fig,'ionflow')
                else
                    plot_ionflow(V_i,absV,T_i,date,expval,ICCD,pathname,factor,mpoints,false,save_fig,'ionflow')
                end
            end
        end
    end
else
    error('begin_cal must <= %d.', exp_log(end_row,4))
end
