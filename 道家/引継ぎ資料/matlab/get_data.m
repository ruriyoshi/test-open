function [data,TF_data] = get_data(date,shot,TF_shot,smoothing,TF_offset,movemean,t_start,t_end)

% folder directory
%folder_directory = 'C:/Users/take_/OneDrive/デスクトップ/program/fluctuation/data/';
%folder_directory = 'C:/Users/take_/OneDrive/デスクトップ/program/モード計測/Data/fluctuation/';
%folder_directory = 'D:\pub\mnt\old-koala\experiment\results\fluctuation\data/';
folder_directory = 'E:/experiment/results/fluctuation/data/';  %old coala
[data_dir,TF_dir] = directory_generation(date,shot,TF_shot);

% coeffiocient
N = 13;
S = 1e-3 * 5e-3;
R = 1e7;
C = 1e-10;


% acquire dB/dt
raw_data = readmatrix(data_dir,'TrimNonNumeric',true);
raw_TF_data = readmatrix(TF_dir,'TrimNonNumeric',true);
data_array = raw_data(4:end,2);
TF_data_array = raw_TF_data(4:end,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%↓
porality = [1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,-1,-1,1,1,1,1,1,1,1,-1,-1];


% create matrix data(time,CH)
% CH... 1 to 24
t = t_start:t_end;
data = zeros(length(t),24);
TF_data = zeros(length(t),24);

for i = 1:24
    data(:,i) = data_array(length(t)*(i-1)+1:i*length(t));
    TF_data(:,i) = TF_data_array(length(t)*(i-1)+1:i*length(t));
end

% transpose and adjust porality
data = data';
data = data.*(porality)';
TF_data = TF_data';
TF_data = TF_data.*(porality)';

% ## ged rid of TF offset ##
if TF_offset
    data = data - TF_data;
end


%calibration
data = data/(N*S);

if smoothing
    for i = 1:24
        data(i,:) = smoothdata(data(i,:),'rloess',15);
    end
end

if movemean
    data = movmean(data,5,2);
end

% assign value to dead channel\
 data(7,:) = 0.5*data(5,:) + 0.5*data(9,:);
 %data(6,:) = 0.5*data(7,:) + 0.5*data(5,:);
 data(:,[13 14]) = data(:,[14 13]);
 %data(18,:) = data(7,:);
 data(17,:) = data(19,:);
 data(19,:) = 0.5*data(18,:) + 0.5*data(20,:);

function [data_dir,TF_dir] = directory_generation(date,shot,TF_shot)
    if shot < 10
        data_dir = [folder_directory,num2str(date),'/',num2str(date),'00',num2str(shot),'.csv'];
    elseif shot < 100
        data_dir = [folder_directory,num2str(date),'/',num2str(date),'0',num2str(shot),'.csv'];
    elseif shot < 1000
        data_dir = [folder_directory,num2str(date),'/',num2str(date),num2str(shot),'.csv'];
    else
        disp('More than 999 shots! You need some rest!!!')
    end
        
    if TF_shot < 10
        TF_dir = [folder_directory,num2str(date),'/',num2str(date),'00',num2str(TF_shot),'.csv'];
    elseif TF_shot < 100
        TF_dir = [folder_directory,num2str(date),'/',num2str(date),'0',num2str(TF_shot),'.csv'];
    elseif TF_shot < 1000
        TF_dir = [folder_directory,num2str(date),'/',num2str(date),num2str(TF_shot),'.csv'];
    else
        disp('More than 999 shots! You need some rest!!!')
    end
end

clearvars -except data TF_data

end