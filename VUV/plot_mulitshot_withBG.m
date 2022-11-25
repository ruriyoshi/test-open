function [f1,f2] = plot_mulitshot_withBG(date,BG,shot,ch,pass)
% You can plot some shot data of oscilloscope that you need using the date and the shot.
% This program is focus on the case the data has some background data.
% You need 'make_data_dir.m' and 'get_one_ch.m' and 'lowpass_fft.m' file to use this function. 
% input 
%     date: the expriment day
%     BG  : arrangement of the background data like this, BG = [1 2 3 4];
%           If you have only one background(the background data is same all
%           time), you can set BG only one shot like this, BG = 1;
%     shot: arrangement of the main shot data like this, shot = [5 6 7 8];
%     ch  : oscilloscope channel you need
%     pass: If you need remove high frequency noise, you can set pass as
%           the cut off frequency like this, pass = 10000; if you don't,
%           you should set it as false. (pass = false;)
% output
%     f1  : Figure object. All data(BG and main shot data) is plotted on this figure.
%     f2  : Figure object. main shot data minus the background data is
%           plotted on this figure.
    data = [];
    data_BG = [];
    for i = BG
        data1 = get_one_ch(date,i,ch);
        if pass ~= false
            [~, data1] = lowpass_fft(data1,pass);
        end
        data_BG = [data_BG data1];
    end
    for i = shot
        data1 = get_one_ch(date,i,ch);
        if pass ~= false
            [~, data1] = lowpass_fft(data1,pass);
        end
        data = [data data1];
    end
    
    f1 = figure;
    hold on
    for i = 1:length(data_BG(1,:))
        txt = ['BG,shot=',num2str(BG(i))];
        plot(data_BG(1:end,i),'DisplayName',txt)
    end
    for i = 1:length(data(1,:))
        txt = ['plasma,shot=',num2str(shot(i))];
        plot(data(1:end,i),'DisplayName',txt)
    end
    legend show
    grid on
    xlim([3953 5763])
    ylim([-0.163 0.088])
    
    f2 = figure;
    hold on
    for i = 1:length(data(1,:))
        txt = ['shot=',num2str(shot(i))];
        if length(BG) == 1
            plot(data(:,i)-data_BG,'DisplayName',txt);
        else
            plot(data(:,i)-data_BG(:,i),'DisplayName',txt);
        end
    end
    legend show
    grid on
    xlim([3953 5763])
    ylim([-0.163 0.088])
end