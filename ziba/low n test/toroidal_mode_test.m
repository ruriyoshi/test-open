clear all
% *************** generate test signal ***************
t = 1:200;
%x = 0:1/8:7/8; % uniform probe placement
x = [0,4/24,7/24,9/24,12/24,15/24,17/24,20/24]; % non-uniform probe placement
% amplitude modulation
n0 = 3;
n1 = 1+0.8*sin(2*pi*t/t(end));
n2 = exp(t/t(end));
n3 = 1+0.6*sin(pi/2+2*pi*t/t(end));
% frequency modulation
f1 = 0.03+0.08*cos(4*pi*t/t(end));
f2 = 0.10+0.08*sin(10*pi*t/t(end));
f3 = 0.18+0.05*sin(6*pi*t/t(end));
% generate signal
test_signal = zeros(length(t),length(x));
for i = 1:length(t)
    test_signal(i,:) = n0 + ...
                  n1(i)*sin(1*2*pi*x+sum(f1(1:i))*2*pi) + ...
                  n2(i)*sin(2*2*pi*x+sum(f2(1:i))*2*pi) + ...
                  n3(i)*sin(3*2*pi*x+sum(f3(1:i))*2*pi) + ...
                  randn(1,length(x))*0.05; % added random noise
end

% *************** calculate mode amplitude and frequency ***************
% signal format: n_time x n_channel (example: 1000 x 8 for real experiment)
% can subsitude test_signal with real signal in real experiment
[n_Amp,n_Ph] = toroidal_mode(t,x,test_signal); % get amplitude and phase

% calculate frequency from phase (smoothed)
n1_Omega = calculate_omega(smooth(cumulative_phase(n_Ph(1,:)),0.1,'loess'));
n2_Omega = calculate_omega(smooth(cumulative_phase(n_Ph(2,:)),0.1,'loess'));
n3_Omega = calculate_omega(smooth(cumulative_phase(n_Ph(3,:)),0.1,'loess'));

% *************** plotting ***************
figure('Position', [0 0 800 800],'visible','on');
subplot(7,4,1:4)
hold on
plot(t,smooth(n_Amp(2,:)./n_Amp(1,:),0.1,'loess'),'k');
plot(t,smooth(n_Amp(3,:)./n_Amp(1,:),0.1,'loess'),'b');
plot(t,smooth(n_Amp(4,:)./n_Amp(1,:),0.1,'loess'),'r');
scatter(t,n_Amp(2,:)./n_Amp(1,:),6,'k');
scatter(t,n_Amp(3,:)./n_Amp(1,:),6,'b');
scatter(t,n_Amp(4,:)./n_Amp(1,:),6,'r');
legend('n=1','n=2','n=3');
ylim([0 1])
title('calculated mode amplitude')
hold off

subplot(7,4,5:8)
hold on
plot(t,n1/n0,'k');plot(t,n2/n0,'b');plot(t,n3/n0,'r');
legend('n=1','n=2','n=3');
title('theoretical mode amplitude')
hold off

subplot(7,4,9:12)
hold on
plot(t,n_Ph(1,:),'k');plot(t,n_Ph(2,:),'b');plot(t,n_Ph(3,:),'r');
ylim([-pi pi]);
legend('n=1','n=2','n=3');
title('calculated phase')
hold off

subplot(7,4,13:16)
hold on
plot(t,rem(total_phase(f1)+pi/2,2*pi)-pi,'k');
plot(t,rem(total_phase(f2)+pi/2,2*pi)-pi,'b');
plot(t,rem(total_phase(f3)+pi/2,2*pi)-pi,'r');
ylim([-pi pi]);
legend('n=1','n=2','n=3');
title('theoretical phase')
hold off

subplot(7,4,17:20)
hold on
plot(t,n1_Omega/(2*pi),'k');
plot(t,n2_Omega/(2*pi),'b');
plot(t,n3_Omega/(2*pi),'r');
scatter(t,calculate_omega(cumulative_phase(n_Ph(1,:)))/(2*pi),6,'k');
scatter(t,calculate_omega(cumulative_phase(n_Ph(2,:)))/(2*pi),6,'b');
scatter(t,calculate_omega(cumulative_phase(n_Ph(3,:)))/(2*pi),6,'r');
legend('n=1','n=2','n=3');
ylim([-0.1 0.3]);
title('calculated frequency')
hold off

subplot(7,4,21:24)
hold on
plot(t,f1,'k');plot(t,f2,'b');plot(t,f3,'r');
legend('n=1','n=2','n=3');
ylim([-0.1 0.3]);
title('theoretical frequency')
hold off

subplot(7,4,25:28)
hold on
for i = 1:8
    plot(t,test_signal(:,i));
    title('test signal')
end
hold off

%Animate_toroidal_mode(t,n_Amp,n_Ph);


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