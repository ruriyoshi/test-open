function [f] = plot_mulitshot(date,shot,ch)
% You can plot some shot data of oscilloscope that you need using the date and the shot.
% You need 'make_data_dir.m' and 'get_one_ch.m' file to use this function. 
% input 
%     date: the expriment day
%     shot: arrangement of the shot like this, [1 2 3 4]
%     ch  : oscilloscope channel you need
% output
%     f   : figure object
    f = figure;
    hold on
    for i = shot
        data = get_one_ch(date,i,ch);
        txt = ['shot=',num2str(i)];
        plot(data,'DisplayName',txt);  
    end
    legend show
    hold off
    xlim([4500 4800]);
    ylim([-0.06 0.02]);
end