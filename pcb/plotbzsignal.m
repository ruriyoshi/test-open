function plotbzsignal(y_upper_lim, col2, col1, t_end, p_ch, y_lower_lim, t_start, bz, ok, r_ch, r)
figure('Position', [10 10 1200 900])
for i=1:r
    for j=1:col1
        subplot(r,col1,(i-1)*col1+j)
        if ok(r_ch*(i-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bz(t_start:t_end,r_ch*(i-1)+j))
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bz(t_start:t_end,r_ch*(i-1)+j),'r:')
        end   
        title(num2str(p_ch(i,j)));
        %xlim([t_start t_end]);
        xticks([t_start t_end]);
        ylim([y_lower_lim y_upper_lim]);
    end
end

figure('Position', [10 10 1200 900])
for i=1:r
    for j=col1+1:col1+col2
        subplot(r,col2,(i-1)*col2+j-col1)
        if ok(r_ch*(i-1)+j)==1 
            plot(t_start:t_end,bz(t_start:t_end,r_ch*(i-1)+j))
        else 
            plot(t_start:t_end,bz(t_start:t_end,r_ch*(i-1)+j),'r:')
        end   
        title(num2str(p_ch(i,j)));
        %xlim([t_start t_end]);
        xticks([t_start t_end]);
        ylim([y_lower_lim y_upper_lim]);
    end
end
end



