function [V_i,absV,T_i] = cal_ionflow(date,ICCD,mpoints,pathname,show_offset,plot_fit,save_fit,save_flow)
%������/ICCD�ϐ�/�v���_�z��/pathname/offset��\��/�K�E�X�t�B�b�e�B���O���v���b�g/������ۑ�

%------�t�B�b�e�B���O�����p�yinput�z-------
w_CH = 4;%�yinput�z�`�����l������(Y����)�������킹��
l_mov = 7;%�yinput�z�g�������ړ����ϒ���
th_ratio = 0.8;%�yinput�z�t�B�b�e�B���O��臒l
l_L1 = 41;%�yinput�z�g�����̐؂��蒷��
d_L1 = 25;%�yinput�z�g�������ʒu���꒲��
d_CH = 2;%�yinput�zCH�����ʒu���꒲��

%�����萔
Vc = 299792.458;%����(km/s)
Angle = 30;%���ˊp(�x)
switch ICCD.line
    case 'Ar'%�A���S���̎�
        A = 40;%���q��
        lambda0 = 480.602;%�g�p�X�y�N�g��(nm)
        lambda1 = 480.7019;%�Z�������v�X�y�N�g��(nm)
        lambda2 = 479.2619;%�Z�������v�X�y�N�g��(nm)
    case 'H'%���f�̎�
        A = 1;%���q��
        lambda0 = 486.135;%�g�p�X�y�N�g��(nm)
        warning('Sorry, not ready for H experiment.')%ICCD.line�̓��̓G���[
        return;
    otherwise
        warning('Input error in ICCD.line.')%ICCD.line�̓��̓G���[
        return;
end
center = load_calibration(date,ICCD);

dir1 = [pathname.IDSP '/' num2str(date)];%�f�B���N�g��1
if mpoints.n_z == 1
    filename1 = [dir1 '/shot' num2str(ICCD.shot) '_' num2str(ICCD.trg) 'us_w=' num2str(ICCD.exp_w) '_gain=' num2str(ICCD.gain) '.asc'];%ICCD�t�@�C����
    if not(exist(filename1,"file"))
        warning([filename1,' does not exist.']);
        V_i = char.empty;
        absV = char.empty;
        T_i = char.empty;
        return
    end
    %data1����X�y�N�g���Q���擾
    data1 = importdata(filename1);
end

%----------------------------------------------------------------------------
% centerX = repmat(center(:,3),1,mpoints.n_z);%�`�����l���Ή����S����X���W
switch mpoints.n_z
    case 1
        centerY = center(:,2);%�`�����l���Ή����SY���W
    case 2
        warning('Sorry, not ready for n_z = 2.')%ICCD.line�̓��̓G���[
        return
        centerY = repmat(center(:,2),1,mpoints.n_z);%�`�����l���Ή����SY���W
end

X1 = data1(:,1);%X1(�s�N�Z��)�����`
[l_X1,~]=size(data1);%X1���̒������擾
L1 = zeros(l_X1,mpoints.n_CH);%L1(�g��)�����`
L1_shaped = zeros(l_L1,mpoints.n_CH);%�؂������g����
px2nm = zeros(mpoints.n_CH,1);%nm/pixel
switch ICCD.line
    case 'Ar'%�A���S���̎�
        %     lambda = [lambda1 lambda2];
        for i = 1:mpoints.n_CH
            %         pixel = [center(i,3) center(i,4)];
            %         p = polyfit(pixel,lambda,1);
            %         px2nm(i,1) = p(1);
            %         L1(:,i) = polyval(p,X1);
            px2nm(i,1) = center(i,4);
            L1(:,i) = lambda1 - px2nm(i,1)*(X1(:,1)-center(i,3));
            L1_shaped(:,i) = L1(round(center(i,3))-(l_L1-1)/2+d_L1:round(center(i,3))+(l_L1-1)/2+d_L1,i);
        end
    case 'H'%���f�̎�
        for i = 1:mpoints.n_CH
            px2nm(i,1) = 0.00536;
            L1(:,i) = px2nm(i,1)*(X1-center(i,3))+lambda0;
        end
end
spectrum1 = zeros(l_L1,mpoints.n_CH);%data1�̕������ʂ�����
for i = 1:mpoints.n_CH
    spectrum1(:,i) = ...
        sum(data1(round(center(i,3))-(l_L1-1)/2+d_L1:round(center(i,3))+(l_L1-1)/2+d_L1,round(centerY(i,1)+d_CH-w_CH):round(centerY(i,1)+d_CH+w_CH)),2);
end

% %data2����X�y�N�g���Q���擾
% if mpoints.n_z > 1
%     data2 = importdata(filename2);
%     X2 = data2(:,1);%X2(�s�N�Z��)�����`
%     [LofX2,~]=size(data2);%X2���̒������擾
%     spectrum2=zeros(LofX2,mpoints.n_CH);%data2�̕������ʂ�����
%     for i = 1:mpoints.n_CH
%         spectrum2(:,i) = ...
%             sum(data2(:,round(centerY(i,2)-width):round(centerY(i,2)+width)),2);
%     end
% end

%�g�������̈ړ����ς��痣�ꂽ�O��l���ړ����ςɕύX(�m�C�Y����)
M1 = movmean(spectrum1,l_mov);
% M1 = (M1*l - spectrum1) / (l-1);
% if mpoints.n_z == 2
%     M2 = movmean(spectrum2,mov_l);
%     % M2 = (M2*l - spectrum2) / (l-1);
% end
for i = 1:mpoints.n_CH
    for j = 1:l_L1
        if spectrum1(j,i) > M1(j,i)*1.1
            spectrum1(j,i) = M1(j,i);
        elseif spectrum1(j,i) < M1(j,i)*0.9
            spectrum1(j,i) = M1(j,i);
        end
        % if mpoints.n_z == 2
        %     if spectrum2(j,i) > M2(j,i)*1.1
        %         spectrum2(j,i) = M2(j,i);
        %     elseif spectrum2(j,i) < M2(j,i)*0.9
        %         spectrum2(j,i) = M2(j,i);
        %     end
        % end
    end
end

%�v�Z���ʎ擾�p�̔z���p��
amp = zeros(mpoints.n_CH,mpoints.n_z);%�U��1���data1�A2���data2
shift = zeros(mpoints.n_CH,mpoints.n_z);%���S1���data1�A2���data2
sigma = zeros(mpoints.n_CH,mpoints.n_z);%�V�O�}1���data1�A2���data2
V_i = zeros(mpoints.n_r,mpoints.n_z*2);%1���data1�EVz(km/s)�A2���data1�EVr(km/s)�A3���data2�EVz(km/s)�A4���data2�EVr(km/s)
absV = zeros(mpoints.n_r,mpoints.n_z);%1���data1�EV(km/s)�A2���data2�EV(km/s)
T_CH = zeros(mpoints.n_CH,mpoints.n_z);%1���data1�E���x(eV)�A2���data1�E���x(eV)
T_i = zeros(mpoints.n_r,mpoints.n_z);%1���data1�E���x(eV)�A2���data1�E���ω��x(eV)
offset = zeros(mpoints.n_r,mpoints.n_z);
error_CH = zeros(mpoints.n_CH,1);

for k = 1:mpoints.n_r
    %0�x�y�A�X�y�N�g������I�t�Z�b�g�����o
    for i = 1:2
        i_CH = (k-1)*4+i;%CH�ԍ�
        S1 = [L1_shaped(:,i_CH) spectrum1(:,i_CH)]; %[�g��,���x]
        s1 = size(S1);%S1��[�s��,��]
        MAX1 = max(spectrum1(:,i_CH)); %�X�y�N�g���̍ő�l
        j = 1;
        while j < s1(1)+1 %SN�̈����f�[�^������
            if S1(j,2) < MAX1*th_ratio
                S1(j,:) = [];
            else
                j = j+1;
            end
            s1 = size(S1);
        end
        try
            f = fit(S1(:,1),S1(:,2),'gauss1');
        catch ME
            warning(['Fitting failed in CH ',num2str(i_CH),'.']);
            error_CH(i_CH,1) = 1;
            break
        end
        coef=coeffvalues(f);
        shift(i_CH,1) = coef(2)-lambda0;
    end
    offset(k,1) = (shift((k-1)*4+1,1) + shift((k-1)*4+2,1))/2;%�Ό��������瓾��ꂽ�I�t�Z�b�g[nm]
end

%�ُ�l������offset����}���Ēu������
ori_offset = offset;%�ύX�O��ۊ�
buf = offset(:,1);%offset���R�s�[
ind_0 = find(buf == 0);%0�����o
ind_non0 = find(buf ~= 0);%non0�����o
buf(ind_0) = [];%0���폜
buf = filloutliers(buf,"linear");%�O��l����^���}�Œu������
offset(ind_non0,1) = buf;%offset���X�V
offset(ind_0,1) = mean(offset(ind_non0,1));%offset���X�V
for k = 1:mpoints.n_r
    if offset(k,1) ~= ori_offset(k,1)
        disp(['Offset in Row ',num2str(k),' is replaced from ',num2str(ori_offset(k,1)),' to ',num2str(offset(k,1))'.'])
    end
end

%data1�̃K�E�X�t�B�b�e�B���O
if plot_fit
    figure('Position',[300 50 1000 1000],'visible','on')
    sgtitle(['Fitting data (Horizontal�FView Line, Vertical�FMeasured Position)',newline, ...
        'shot',num2str(ICCD.shot),'-',num2str(ICCD.trg),'us-w=',num2str(ICCD.exp_w),'-gain=',num2str(ICCD.gain),'.asc'])
end
for k = 1:mpoints.n_r
    % if (error_CH((k-1)*4+1,1) == 0) && (error_CH((k-1)*4+2,1) == 0)%�΍R�����̃t�B�b�e�B���O�ɐ���
        if show_offset
            disp(['Offset = ',num2str(offset(k,1)),' in Row ',num2str(k),'.'])
        end
        %�I�t�Z�b�g�������ăK�E�X�t�B�b�e�B���O
        for i = 1:4
            i_CH = (k-1)*4+i;%CH�ԍ�
            S1 = [L1_shaped(:,i_CH)-offset(k,1) spectrum1(:,i_CH)]; %[�g��,���x]
            s1 = size(S1);%S1��[�s��,��]
            MAX1 = max(spectrum1(:,i_CH)); %�X�y�N�g���̍ő�l
            j = 1;
            while j < s1(1)+1 %SN�̈����f�[�^������
                if S1(j,2) < MAX1*0.5
                    S1(j,:) = [];
                else
                    j = j+1;
                end
                s1 = size(S1);
            end
            try
                f = fit(S1(:,1),S1(:,2),'gauss1');
                coef=coeffvalues(f);
                amp(i_CH,1) = coef(1);
                shift(i_CH,1) = coef(2)-lambda0;
                sigma(i_CH,1) = sqrt(coef(3)^2-(center(i_CH,5)*px2nm(i_CH,1))^2);
                T_CH(i_CH,1) = 1.69e8*A*(2*sigma(i_CH,1)*sqrt(log(2))/lambda0)^2;
                if plot_fit
                    subplot(mpoints.n_r,4,i_CH);
                    plot(f,S1(:,1),S1(:,2));
                    % plot(S1(:,1),S1(:,2));
                    xline(lambda0);
                    title([num2str(T_CH(i_CH,1)),' eV'])
                    legend('off')
                    xlabel('Wavelength [nm]')
                    ylabel('Intensity [cnt]')
                end
            catch ME
                warning(['Fitting failed in CH ',num2str(i_CH),'.']);
                break
            end
        end
    % else
    %     warning(['Offset failed in Row ',num2str(k),'.']);
    % end
end
if plot_fit
    if save_fit
        time = round(ICCD.trg+ICCD.exp_w/2);%�v������
        if not(exist([pathname.fig,'/fit/',num2str(date)],'dir'))
            mkdir(sprintf("%s/fit", pathname.fig), sprintf("%s", num2str(date)));
        end
        saveas(gcf,[pathname.fig,'/fit/',num2str(date),'/','shot', num2str(ICCD.shot),'_',num2str(time),'us_fit.png'])
        hold off
        close
    else
        hold off
    end
end

%data1�̗����A���x���v�Z
for i = 1:mpoints.n_r
    set = (i-1)*4;
    Va = -shift(set+4,1)/lambda0*Vc;
    Vb = -shift(set+3,1)/lambda0*Vc;
    V_i(i,1) = (Va-Vb)/(2*cos(Angle*pi/180));%Vz
    V_i(i,2) = -(Va+Vb)/(2*sin(Angle*pi/180));%Vr
    absV(i,1) = sqrt(V_i(i,1)^2 + V_i(i,2)^2);
    %        fprintf('(z, r) = (%.2f, %.2f) cm\n (Vz, Vr) = (%.2f, %.2f) km/s\n'...
    %             ,z(i,1),r(i,1),V(i,1),V(i,2));
    %             fprintf('T_i: %.2f, %.2f, %.2f, %.2f eV\n'...
    %                     ,T(set+1,1),T(set+2,1),T(set+3,1),T(set+4,1));
    %         T_i(i,1) = (T(set+1,1)+T(set+2,1)+T(set+3,1)+T(set+4,1))/4;
    T_i(i,1) = trimmean(T_CH(set+1:set+4,1),50);
    %         fprintf('%.2f\n',T_i(i,1));
end

% %data2�̃K�E�X�t�B�b�e�B���O
% if mpoints.n_z == 2
%     if plot_fit
%         figure('Position',[300 50 1000 1000])
%     end
%     for k = 1:mpoints.n_r
%         %0�x�y�A�X�y�N�g������I�t�Z�b�g�����o
%         for i = 1:2
%             i_CH = (k-1)*4+i;%CH�ԍ�
%             S2 = [L1_shaped(:,i_CH) spectrum2(:,i_CH)]; %[�g��,���x]
%             s2 = size(S1);%S1��[�s��,��]
%             MAX2 = max(spectrum2(:,i_CH)); %�X�y�N�g���̍ő�l
%             j = 1;
%             while j < s2(1)+1 %SN�̈����f�[�^������
%                 if S2(j,2) < MAX2*0.5
%                     S2(j,:) = [];
%                 else
%                     j = j+1;
%                 end
%                 s2 = size(S2);
%             end
%             f = fit(S2(:,1),S2(:,2),'gauss1');
%             coef=coeffvalues(f);
%             amp(i_CH,2) = coef(1);
%             shift(i_CH,2) = coef(2)-lambda0;
%             sigma(i_CH,2) = sqrt(coef(3)^2-(center(i_CH,5)*px2nm(i_CH,1))^2);
%             T_CH(i_CH,2) = 1.69e8*A*(2*sigma(i_CH,2)*sqrt(log(2))/lambda0)^2;
%         end
%         offset(k,2) = (shift((k-1)*4+1,2) + shift((k-1)*4+2,2))/2;%�Ό��������瓾��ꂽ�I�t�Z�b�g[nm]
%         if show_offset
%             disp(offset(k,2))
%         end
%         %�I�t�Z�b�g�������ăK�E�X�t�B�b�e�B���O
%         for i = 1:4
%             i_CH = (k-1)*4+i;%CH�ԍ�
%             S2 = [L1_shaped(:,i_CH) spectrum2(:,i_CH)]; %[�g��,���x]
%             s2 = size(S1);%S1��[�s��,��]
%             MAX2 = max(spectrum2(:,i_CH)); %�X�y�N�g���̍ő�l
%             j = 1;
%             while j < s2(1)+1 %SN�̈����f�[�^������
%                 if S2(j,2) < MAX2*0.5
%                     S2(j,:) = [];
%                 else
%                     j = j+1;
%                 end
%                 s2 = size(S2);
%             end
%             f = fit(S2(:,1),S2(:,2),'gauss1');
%             coef=coeffvalues(f);
%             amp(i_CH,2) = coef(1);
%             shift(i_CH,2) = coef(2)-lambda0;
%             sigma(i_CH,2) = sqrt(coef(3)^2-(center(i_CH,5)*px2nm(i_CH,1))^2);
%             T_CH(i_CH,2) = 1.69e8*A*(2*sigma(i_CH,2)*sqrt(log(2))/lambda0)^2;
%             if plot_fit
%                 subplot(mpoints.n_r,4,i_CH);
%                 plot(f,S2(:,1),S2(:,2));
%                 xline(lambda0);
%                 title([num2str(T_CH(i_CH,2)),' eV'])
%                 legend('off')
%                 xlabel('Wavelength [nm]')
%                 ylabel('Intensity [cnt]')
%             end
%         end
%     end
%     if plot_fit
%         sgtitle('Fitting data1 (Horizontal�FChannel Vertical�FPosition)')
%     end
%
%     %data2�̗����A���x���v�Z
%     for i = 1:mpoints.n_r
%         set = (i-1)*4;
%         Va = -shift(set+4,2)/lambda0*Vc;
%         Vb = -shift(set+3,2)/lambda0*Vc;
%         V_i(i,3) = (Va-Vb)/(2*cos(Angle*pi/180));%Vz
%         V_i(i,4) = -(Va+Vb)/(2*sin(Angle*pi/180));%Vr
%         absV(i,2) = sqrt(V_i(i,3)^2 + V_i(i,4)^2);
%         %        fprintf('(z, r) = (%.2f, %.2f) cm\n (Vz, Vr) = (%.2f, %.2f) km/s\n'...
%         %             ,z(i,2),r(i,2),V(i,3),V(i,4));
%         %             fprintf('T_i: %.2f, %.2f, %.2f, %.2f eV\n'...
%         %                     ,T(set+1,2),T(set+2,2),T(set+3,2),T(set+4,2));
%         T_i1=(T_CH(set+1,2)+T_CH(set+2,2)+T_CH(set+3,2)+T_CH(set+4,2))/4;
%         fprintf('%.2f\n',T_i1);
%     end
% end

%�����f�[�^��ۑ�
if save_flow
    if not(exist([pathname.flowdata,'/',num2str(date)],'dir'))
        mkdir(sprintf("%s", pathname.flowdata), sprintf("%s", num2str(date)));
    end
    save([pathname.flowdata,'/',num2str(date),'/shot',num2str(ICCD.shot),'_',num2str(ICCD.trg),'us_w=',num2str(ICCD.exp_w),'_gain=',num2str(ICCD.gain),'.mat'],'V_i','absV','T_i')
end