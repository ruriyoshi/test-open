function [dBdt] = contour_fluctuation(data,t_start,t_end,shot,set_caxis)

% set axis
ch_position = -172.5:15:172.5;
t = linspace(t_start,t_end,t_end*10-t_start*10+1);
t_ticks = linspace(t_start,t_end,t_end-t_start+1);

ch_tick = zeros(1,24);
for i=1:24
    ch_tick(i) = i;
end
str_ch_position = string(ch_position);
y_tick = "ch" + string(ch_tick) + "  " + str_ch_position;


% extract data
dBdt = data(:,t_start*10:t_end*10);
%dBdt = data;

% plot 
figure('Position', [10 10 800 400],'name',['shot', num2str(shot)]);
contourf(t,ch_position,dBdt,10,'LineStyle','none');
title('fluctuations');
xlim([t_start t_end]);
%xticks(t_ticks);
xlabel('time[us]','Fontsize',18,'FontWeight','bold');
ylabel('toroidal angle [rad]','Fontsize',18,'FontWeight','bold');
g = gca;
g.XAxis.FontSize = 12;
g.YAxis.FontSize = 12;
g.ZAxis.FontSize = 12;
yticks([-172.5, 0, 172.5]);
yticklabels({'0','π/4','π/2'});
c = colorbar;
if set_caxis
    caxis([-1000 1000]);
end
colormap turbo;
%colorbar('Fontsize',11,'FontWeight','bold')
c.Label.String = 'dB/dt [T/s]';
c.Label.FontSize = 15;
c.Label.FontWeight = 'bold';
%caxis([-0.15 0.15]);


end