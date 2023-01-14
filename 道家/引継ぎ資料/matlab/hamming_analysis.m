
t = linspace(430,520,901);
l = length(t);

%figure('name','raw data');
%plot(t,data(12,4300:5200));

data_hamming = zeros(size(data)); % data_hamming(24,10001)   0
hamming_window = zeros(1,10001);  % hamming_window(1,10001)  0
window = (hamming(l))';           % window(1,901)  value
hamming_window(1,4300:5200) = window; % hamming_window(1,901) 
data_hamming(:,4300:5200) = data(:,4300:5200).*hamming_window(1,4300:5200);

% [power_rec] = power_spectrum_fluctuation(data,shot,4300,5200,12);
% [power_hamming] = power_spectrum_fluctuation(data_hamming,shot,4300,5200,12);
% 
% 
% figure;
% plot(data_hamming(12,4300:5200));
% legend('hamming data')
% 
% figure;
% plot(window);
% legend('hamming window');
% 
