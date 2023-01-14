clear all

date = 220103;
TF_shot = 2 ;
TF_offset = true;
shot = 39 ;

% acquire data...get_data(date,sh51ot,TF_shot,smoothing,TF_offset,movemean,t_start,t_end)
[data,TF_data] = get_data(date,shot,TF_shot,false,TF_offset,true,0,10000);


% plot dB/dt...contourf_fluctuation(data,t_start,t_end,set_caxis)
%[dBdt] = contour_fluctuation(data,430,520,shot,true);

% ##plot signal(data, move_mean) ##
%plot_signal(data,true);

% ##toroidal mode##
%parameters(data,t_start,t_end,shot,show_mode);
%[amp,toroidal_mode] = high_n_toroidal_mode(data,460,490,shot,true);


% ## plot power spectrum ##
%parameters(data,shot,t_start,t_end,ch_number)
%[power] = power_spectrum_fluctuation(data,shot,4300,5200,21);
%[spectrum,ifft_spectrum] = old_delta_B(data,shot,4300,6000,1);
%[delta_B,frequency] = new_delta_B(data,shot,430,520,1e5,5);


% ## animation ##
%anime_fluctuation(data);

