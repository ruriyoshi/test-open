%function [delta_B,frequency] = new_delta_B(data,shot,t_start,t_end,filter,window_size)
t_start = 450;
t_end = 520;
window_size = 20;
st = t_start*10-1;
filter = 2e5;

window = window_size*10;
hann_window = (hann(window))';
fs = 1e7;
x = linspace(t_start,t_end,t_end*10-t_start*10);
delta_B = zeros(24,length(x));
B = zeros(24,length(x));
delta = zeros(24,length(x));

ch_tick = zeros(1,24);
for i=1:24
    ch_tick(i) = i;
end
ch_position = -172.5:15:172.5;
str_ch_position = string(ch_position);
y_tick = "ch" + string(ch_tick) + "  " + str_ch_position;


% hht analysis
for t = t_start*10:t_end*10
    
    t_window_st = t - fix(window/2);
    t_window_end = t + fix(window/2) -1;
    
    hann_data = data(:,t_window_st:t_window_end).*hann_window;
    fft_data = fft(hann_data(:,:),[],2);
    l = length(fft_data);
    fft_data = fft_data(:,1:l/2+1);
    fft_data(:,2:end-1) = 2*fft_data(:,2:end-1);
    frequency = fs*(0:(l/2))/l;
    frequency(1,1) = 1;
    omega = 2*pi*frequency;
    
    
%     filtering
      filter_index = frequency >= filter;
      filter_data = fft_data.*filter_index;
     

    nolm_delta = abs(filter_data./omega);
    nolm_B = abs(fft_data./omega);
    %nolm_B(:,1) = 1;

    sum_delta = sum(nolm_delta,2);
    sum_B = sum(nolm_B,2);
    
    B(:,t-st) = sum_B;
    delta(:,t-st) = sum_delta;
    
    delta_B(:,t-st) = sum_delta./sum_B;
end
    

figure('Position', [10 10 800 400],'name',['shot', num2str(shot)]);
contourf(x,ch_position,delta(:,1:end-1),15,'LineStyle','none');
title('fluctuations');
%xlim([t_start t_end]);
%xticks(t_ticks);
xlabel('time[us]','Fontsize',18,'FontWeight','bold');
ylabel('toroidal position[mm]','Fontsize',18,'FontWeight','bold');
yticks(ch_position);
yticklabels(y_tick);
c = colorbar;
% if set_caxis
%     caxis([-700 800]);
% end
colormap turbo;
%end
