function return_data = get_one_ch(date,shot,ch_num)
% You can get one channel data of oscilloscope that you need using the date and the shot.
% You need make_data_dir.m file to use this function. 
    data_dir = make_data_dir(date,shot);
    data = dlmread(data_dir,'\t',1,1);
    return_data = data(:,ch_num+1);
end