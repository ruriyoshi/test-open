function [Et_x_point] = plot_Et_in_time_at_x_point(B_z,r_probe,t_start,t_end,shot)

 [x_point_index] = x_point(B_z,r_probe);
 [Et] =  get_Et(B_z,r_probe);
 
 Et(:,:,1024) = Et(:,:,1023);
 Et_x_point = (Et(x_point_index))';
 
 t = linspace(t_start,t_end,t_end-t_start+1);
 
 
 figure('name',['shot', num2str(shot)]);
 
 plot(t,-Et_x_point(1,t_start:t_end));
 ylim([-50 250]);
% xlim([455 490]);
 g = gca;
 g.XAxis.FontSize = 12;
 g.YAxis.FontSize = 12;
 title('-Et at x point');
 xlabel('time [us]','Fontsize',18,'FontWeight','bold');
 ylabel('-Et [V/m]','Fontsize',18,'FontWeight','bold');
 
end
 
 