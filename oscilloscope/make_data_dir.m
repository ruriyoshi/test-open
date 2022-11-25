function data_dir = make_data_dir(date,shot)
% You can get the directory path you need using the date and the shot.
% You must change the part '/Users/keisukemiki/koala_mnt' to your directory
%    full path.
    if shot < 10
        data_dir = ['/Users/keisukemiki/koala_mnt/experiment/results/ts-3u/',num2str(date),'/',num2str(date),'00',num2str(shot),'.rgw'];
    elseif shot < 100
        data_dir = ['/Users/keisukemiki/koala_mnt/experiment/results/ts-3u/',num2str(date),'/',num2str(date),'0',num2str(shot),'.rgw'];
    elseif shot < 1000
        data_dir = ['/Users/keisukemiki/koala_mnt/experiment/results/ts-3u/',num2str(date),'/',num2str(date),num2str(shot),'.rgw'];
    else
        disp('More than 999 shots! You need some rest!!!')
    end
end