function [power] = power_spectrum_fluctuation(data,shot,t_start,t_end,ch_num)

fs = 1e7;
%data = data(:,t_start:t_end);
t = linspace(t_start,t_end,t_end-t_start);
l = length(t);

spectrum = fft(data,[],2);
p = abs(spectrum).^2/l;
power = p(:,1:l/2+1);
power(:,2:end-1) = 2*power(:,2:end-1);
frequency = fs*(0:(l/2))/l;


figure('name',['shot', num2str(shot)]);
plot(frequency,power(ch_num,:));
%semilogx(frequency,power(ch_num,:));
xlabel('frequency [Hz]');
ylabel('power');
title('power spectrum');
legend("ch" + num2str(ch_num));

end