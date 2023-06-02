function [Bt, Bz] = angle_calibration(t,test_Bt_raw, test_Bz_raw)

probe_num = 1;


% signal format: n_time x n_channel (example: 1000 x 8 for real experiment)
% can subsitude test_signal with real signal in real experiment

%initialization
    %alpha = zeros(1, probe_num);
    theta = zeros(1, probe_num);
    %one_list = ones(1, probe_num);
    B = zeros(2,length(t), probe_num);
    

for i = 1:length(t)
    
    %alpha(i,:) = test_Bt_raw(i,:)/test_Bz_raw(i,:);% 1*channnel_num
    %theta = arctan(1/alpha)

    theta(i,:) = atan2(test_Bz_raw(i,:),test_Bt_raw(i,:)); % 1*channel_num
  

    for j = 1:probe_num
        mat = [cos(theta(i,j)) sin(theta(i,j)); -sin(theta(i,j)) cos(theta(i,j))];
        B(:,i,j) = mat*[test_Bt_raw(i,j); test_Bz_raw(i,j)];
        Bt(i,j) = B(1,i,j);
        Bz(i,j) = B(2,i,j);
    end

end

end

