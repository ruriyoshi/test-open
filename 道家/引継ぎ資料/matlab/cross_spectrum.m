function cross_spectrum(data,t_start,t_end,ch_a,ch2_b)

cross_data = data(:,t_start:t_end);

fs = 1e7;
t = linspace(t_start,t_end,(t_end-t_start)*10+1);
l = length(t);
end