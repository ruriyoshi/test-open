tf = 2.0:0.2:3.4;
before = [478.3, 480.2, 477.2, 476.8, 477, 478.4, 476.3, 476.5, 476.9, 476.8, 476.9, 475.4, 475.1, 474.5, 473.6, 474.3];
after = [474.5, 477.2, 474.3, 473.5, 474.5, 475.7, 474.1, 474.4, 475.4, 474.8, 475.8, 473.8, 474.2, 473.7, 472.8, 473.3];


div = before - after;
div = div * 1e-6;
div_val = 0.345*1e-3./div;
value = zeros(1,8);
error = zeros(1,8);
for i=1:8
    value(i) = (div_val(2*i-1)+div_val(2*i))/2;
    error(i) = abs(value(i)-div_val(2*i)); 
end

B = linspace(0.1,0.17,8);
alfven = 2.2e13*B*sqrt(1/1e20);

figure;
errorbar(tf,value,error);
hold on;
plot(tf,alfven);
g = gca;
g.XAxis.FontSize = 12;
g.YAxis.FontSize = 12;
title('toroidal beam');
xlabel('TF [kV]','Fontsize',18,'FontWeight','bold');
ylabel('velosity [km/s]','Fontsize',18,'FontWeight','bold');
legend('experimental velosity','Alfven speed','location','northwest');
hold off;

