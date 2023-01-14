function plot_signal(data,move_mean)
t_start = 3500;
t_end = 6000;
t_number = t_end - t_start + 1;
x = linspace(t_start*0.1,t_end*0.1,t_number);

%## smoothing ##
if move_mean
    data = movmean(data,11,2);
end


for i=0:5
    figurename = "oscilloscope " +  num2str(i+1);
    ch1_name = "ch" + num2str(i*4+1);
    ch2_name = "ch" + num2str(i*4+2);
    ch3_name = "ch" + num2str(i*4+3);
    ch4_name = "ch" + num2str(i*4+4);
    
    
    data((i*4+1),:) = data((i*4+1),:) + 900;
    data((i*4+2),:) = data((i*4+2),:) + 300;
    data((i*4+3),:) = data((i*4+3),:) - 300;
    data((i*4+4),:) = data((i*4+4),:) - 900;

    
    figure('Name',figurename);
    plot(x,data((i*4+1),t_start:t_end));
    ylim([-1500 1500])
    %xlim([400 600]);
    xlabel('time[us]','Fontsize',18,'FontWeight','bold');
    ylabel('signal','Fontsize',18,'FontWeight','bold');
    hold on
    plot(x,data((i*4+2),t_start:t_end));
    plot(x,data((i*4+3),t_start:t_end));
    plot(x,data((i*4+4),t_start:t_end));
    legend(ch1_name,ch2_name,ch3_name,ch4_name);
end


end