function [amp, toroidal_mode] = high_n_toroidal_mode(data,t_start,t_end,shot,show_mode,Et_x_point)

ch_position = -172.5:15:172.5;

t = linspace(t_start,t_end,t_end*10-t_start*10+1);
x = length(ch_position); %24


fft_data = (data)'; %data(10001,24)
y_fft = fft(fft_data,[], 2);

amp = abs(y_fft(:,x/2+1:end));%amp(10001,12)

n = (0:11)*4;

%mode_array = amp.*n;

toroidal_mode = (sum(amp(:,2:end),2))'; %toroidal_mode(1,10001)

figure('name',['shot', num2str(shot)]);
plot(t,toroidal_mode(1,t_start*10:t_end*10));
 g = gca;
 g.XAxis.FontSize = 12;
 g.YAxis.FontSize = 12;
xlabel('time [us]','Fontsize',14,'FontWeight','bold');
ylabel('mode amplitude [T/s]','Fontsize',14,'FontWeight','bold');
ylim([0 3e4]);


% plot each n
if show_mode
    figure('name',['shot', num2str(shot)]);
    subplot(4,1,1);

    plot(t,amp(t_start*10:t_end*10,1));
    hold on;
    plot(t,amp(t_start*10:t_end*10,2));
    plot(t,amp(t_start*10:t_end*10,3));
    plot(t,amp(t_start*10:t_end*10,4));
    ylim([0 3000]);
    g = gca;
    g.XAxis.FontSize = 10;
    g.YAxis.FontSize = 10;
    xlabel('time [us]','Fontsize',11,'FontWeight','bold');
    %ylabel('mode amplitude [T/s]','Fontsize',9,'FontWeight','bold');
    lgd = legend;
    lgd.FontSize = 10;
    lgd.NumColumns = 2;
    legend('n=0','n=4','n=8','n=12');
    hold off;

    subplot(4,1,2);

    plot(t,amp(t_start*10:t_end*10,5));
    hold on;
    plot(t,amp(t_start*10:t_end*10,6));
    plot(t,amp(t_start*10:t_end*10,7));
    plot(t,amp(t_start*10:t_end*10,8));
    ylim([0 3000]);
    g = gca;
    g.XAxis.FontSize = 10;
    g.YAxis.FontSize = 10;
    xlabel('time [us]','Fontsize',11,'FontWeight','bold');
    ylabel('mode amplitude [T/s]','Fontsize',14,'FontWeight','bold');
    lgd = legend;
    lgd.FontSize = 10;
    lgd.NumColumns = 2;
    legend('n=16','n=20','n=24','n=28');
    hold off;

    subplot(4,1,3);

    plot(t,amp(t_start*10:t_end*10,9));
    hold on;
    plot(t,amp(t_start*10:t_end*10,10));
    plot(t,amp(t_start*10:t_end*10,11));
    plot(t,amp(t_start*10:t_end*10,12));
    ylim([0 8000]);
    g = gca;
    g.XAxis.FontSize = 10;
    g.YAxis.FontSize = 10;
    xlabel('time [us]','Fontsize',11,'FontWeight','bold');
    %ylabel('mode amplitude [T/s]','Fontsize',9,'FontWeight','bold');
    lgd = legend;
    lgd.FontSize = 10;
    lgd.NumColumns = 2;
    legend('n=32','n=36','n=40','n=44');
    hold off;
    
    %　論文用
    subplot(4,1,4);
    plot(460:490,-Et_x_point(1,460:490));
    ylim([-50 250]);
    xlim([460 490]);
    g = gca;
    g.XAxis.FontSize = 10;
    g.YAxis.FontSize = 10;
    title('-Et at x point');
    xlabel('time [us]','Fontsize',11,'FontWeight','bold');
    ylabel('-Et [V/m]','Fontsize',14,'FontWeight','bold');
    
    
end


end