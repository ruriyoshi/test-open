function [] = plot_B_z_vs_z(B_z,t,r_probe,z_probe,ch_dist)

figure
hold on
for i = 1:length(r_probe)
    
    plot(z_probe,B_z(i,:,t));
    a = ch_dist(i,:)';
    b = num2str(a);
    c = cellstr(b);
    dx = 0.00; dy = 0.001; % displacement so the text does not overlay the data points
    text(z_probe+dx, B_z(i,:,t)+dy, c);
    xlabel('r (m)','FontSize',18);
    ylabel('Bz(T)','FontSize',18);
    ax = gca;
    ax.FontSize = 18;
end
legend
title(strcat(num2str(t),' us'));
hold off

end