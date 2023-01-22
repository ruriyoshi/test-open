function [] = plot_B_z_vs_r(B_z,t,r_probe,z_probe,ch_dist)

figure
hold on
ZERO = zeros(length(r_probe),1);
plot(r_probe,ZERO,'k');
i_start = 1;
i_end = 28;
for i = i_start:1:i_end
    
    hold on
    B_z_smooth = smooth(r_probe,B_z(:,i,t),'lowess');
    plot(r_probe,B_z_smooth,'k');
    scatter(r_probe,B_z(:,i,t),80,'k','x');
    xlabel('r (m)','FontSize',18);
    ylabel('Bz(T)','FontSize',18);
    ax = gca;
    ax.FontSize = 18;    
    
    %plot(r_probe,B_z(:,i,t));
    a = ch_dist(:,i);
    b = num2str(a);
    c = cellstr(b);
    dx = 0.00; dy = 0.001; % displacement so the text does not overlay the data points
    text(r_probe+dx, B_z(:,i,t)+dy, c);
    
    hold off
end

ylim([-0.1 0.1]);

title(strcat(num2str(t),' us'));
ylabel('Bz (T)','FontSize',18);
xlabel('r (m)','FontSize',18);
hold off
end