function [low_n_signal] = toroidal_mode_offset_new(low_n_data,shot,offset,smoothing,movemean,standarization)

t = linspace(0,1000,10001);


%プロット表示する時間範囲
t_start = 400;
t_end = 600;

porality = [-1, 1, 1, -1, -1, 1, 1, 1];
%probeの挿入角度
not_rad_x = [8.6, 68, 113, 143, 188, 233, 263, 323];
x = not_rad_x*2*pi/360;
ch_num = length(x);
NS_z = 31.8e-3;
NS_t = 35.3e-3;
NS_r = 107e-3;

R = 100e3;
C = 100e-9; 

%時間(s)*チャンネル数(8ch)
low_n_signal = zeros(length(t),ch_num);

%read data from rgw file → low_n_data, low_n_signalにBzの値が入っている
%low_n_dataに全てのチャンネルの

for i=11:18
    low_n_signal(:,i-10) = low_n_data(:,i);
end

%%%%%%%  asign ch port %%%%%%%
%signal(:,1) = (data(109,:))';
% signal(:,2) = (data(96,:))';
% signal(:,5) = (data(99,:))';
% signal(:,8) = (data(98,:))';

% adjust polarity
low_n_signal = low_n_signal.*porality;

% calibration
low_n_signal = (R*C)/(NS_z)*low_n_signal;



% smooth data
if smoothing
    for i=1:8
    low_n_signal(:,i) = smoothdata(low_n_signal(:,i),'Samplepoints',1);
    end
end

% get rid of offset 
if offset
    for i=1:8
     low_n_signal(:,i) = low_n_signal(:,i) - mean(low_n_signal(100:400,i));
    end
end


% ## movemean ##
if movemean
    low_n_signal = movmean(low_n_signal,10,1);
end

% ## get max value ##
if standarization
    max_array = zeros(1,ch_num);
    for i = 1:8
        max_array(1,i) = max(low_n_signal(4200:end,i));
    end
    max_coefficient = max(max(low_n_signal(4200:end,:)));
    low_n_signal = low_n_signal.*(max_coefficient./max_array);
end

% execute fourie transform
[n_Amp, n_Ph] = toroidal_mode(t,x, low_n_signal);

n0_Am = n_Amp(1,:);
n1_Am = n_Amp(2,:);
n2_Am = n_Amp(3,:);
n3_Am = n_Amp(4,:);

n1_Omega = calculate_omega(smooth(cumulative_phase(n_Ph(1,:)),0.1,'loess'));
n2_Omega = calculate_omega(smooth(cumulative_phase(n_Ph(2,:)),0.1,'loess'));
n3_Omega = calculate_omega(smooth(cumulative_phase(n_Ph(3,:)),0.1,'loess'));

%{
%%%%%% older mode analyze code %%%%%%
y_shift = zeros(length(t),length(x));
for i = 1:length(t)
  y_fft = fft(low_n_signal(i,:));
  %fshift = (-n/2:n/2-1);
  y_shift(i,:) = fftshift(y_fft); %周期を0中心にする
end


%  plotting
Am0 = rot90(abs(y_shift(:,5)));
Am1 = rot90(abs(y_shift(:,6)));
Am2 = rot90(abs(y_shift(:,7)));
Am3 = rot90(abs(y_shift(:,8)));

ph1 = angle(rot90(y_shift(:,6)));
ph2 = angle(rot90(y_shift(:,7)));
ph3 = angle(rot90(y_shift(:,8)));

omega1 = calculate_omega(ph1);
omega2 = calculate_omega(ph2);
omega3 = calculate_omega(ph3);
%}


% plotting function

% mode amplitudeのプロット 
figure
ax = gca;
ax.FontSize = 12;
plot(t,n0_Am,'k');
xlim([420 580])
ylim([0 1]);
xlabel('time[μs]','FontSize',11,'FontWeight','bold');
ylabel('Amplitude','FontSize',11,'FontWeight','bold');
hold on
plot(t,n1_Am);
plot(t,n2_Am);
plot(t,n3_Am);
legend('n=0','n=1','n=2','n=3');
%title('amplitude')
hold off


%ch1~ch8までの磁場信号の時間変化を2×4で個々に表示
f = figure('name',['shot', num2str(shot)]);
f.Position(3) = 800;
for i = 1:8
    subplot(2,4,i)
    plot(t,low_n_signal(:,i));
    ax = gca;
    ax.FontSize = 12;
    if i == 1
        ylabel('B[T]','Fontsize',11,'FontWeight','bold');
    end
    xlim([350 600]);
    ylim([-0.2 0.2]);
    title(['ch',num2str(i)]);
end



% subplot(5,4,[1,4])
% hold on
% plot(t,Am1./Am0,'k')
% xlim([t_start t_end])
% ylim([0 2]);
% plot(t,Am2./Am0,'b')
% plot(t,Am3./Am0,'r')
% legend('n=1','n=2','n=3');
% title('amplitude')
% hold off
% 
% 
% subplot(5,4,[5,8])
% hold on
% plot(t,ph1,'k')
% xlim([t_start t_end]);
% plot(t,ph2,'b')
% plot(t,ph3,'r')
% legend('n=1','n=2','n=3');
% title('phase')
% hold off
% 
% subplot(5,4,[9,12])
% hold on
% plot(t,omega1/(2*pi),'k')
% xlim([t_start t_end]);
% ylim([-0.1 0.3]);
% plot(t,omega2/(2*pi),'b')
% plot(t,omega3/(2*pi),'r')
% legend('n=1','n=2','n=3');
% title('omega')
% hold off


% for i = 1:8
%     subplot(5,4,i+12)
%     plot(t,signal(:,i));
%     xlim([t_start t_end]);
%     ylim([-0.1 0.1]);
%     title(num2str(i),'T');
% end

%ch1~ch8までの磁場信号の時間変化重ねて表示
figure('name',['shot', num2str(shot)]);
for i=1:8
    plot(t,low_n_signal(:,i));
    if i == 1
        hold on
    end
    xlabel('time[μs]','FontSize',11,'FontWeight','bold');
    ylabel('B[T]','FontSize',11,'FontWeight','bold');
    ylim([-0.2 0.2]);
    xlim([t_start t_end])
    ax = gca;
    ax.FontSize = 12;
end
legend('CH1','CH2','CH3','CH4','CH5','CH6','CH7','CH8');
hold off


function Omega = calculate_omega(Ph)
    Omega = zeros(size(Ph));
    for j = 1:length(Ph)-1
        Omega(j) = Ph(j+1)-Ph(j);
    end
    Omega(end) = Omega(end-1);
end




function phase = total_phase(freq)
% get phase from frequency
    phase = zeros(size(freq));
    for k = 1:length(phase)
        phase(k) = sum(freq(1:k))*2*pi;
    end
end

function phase = cumulative_phase(Ph)
% get cumulative phase from phase
    phase = zeros(size(Ph));
    phase_reference = 0;
    phase(1) = phase_reference + Ph(1);
    for k = 2:length(phase)
        if Ph(k)-Ph(k-1) > pi
            phase_reference = phase_reference - 2*pi;
        elseif Ph(k)-Ph(k-1) < -pi
            phase_reference = phase_reference + 2*pi;
        end
        phase(k) = phase_reference + Ph(k);
    end
end

%function [] = Animate_toroidal_mode(t,n_Amp,n_Ph)
%     figure
%     xx = linspace(0,1,30)*2*pi;
%     for j = t
%         jj = j - (t(1)-1);
%         subplot(3,1,1)
%         polarplot(xx,n_Amp(2,jj)/n_Amp(1,jj)*sin(1*xx+n_Ph(1,jj))+1);
%         rlim([0 1.5]);
%         subplot(3,1,2)
%         polarplot(xx,n_Amp(3,jj)/n_Amp(1,jj)*sin(2*xx+n_Ph(2,jj))+1);
%         rlim([0 1.5]);
%         subplot(3,1,3)
%         polarplot(xx,n_Amp(4,jj)/n_Amp(1,jj)*sin(3*xx+n_Ph(3,jj))+1);
%         rlim([0 1.5]);
%         
%         pause(0.02);
%     end
% end

end