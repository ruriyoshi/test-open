%function [Bt, Bz] = angle_calibration(Bt_raw, Bz_raw)
clear all

% *************** generate test signal ***************
%direction_num : x, y, z direction → 3

t = 1:200;
T = length(t);
total_channel_num = 24;
probe_num = 8;

direction_num = total_channel_num/probe_num; 


test_Bt_raw = zeros(length(t), probe_num); % 200*8
test_Bz_raw = zeros(length(t), probe_num);


for i = 1:length(t)
   
    offset = 0.05;

    %test_Bt_raw(i,:) = sin(2*pi*i/T); test_Bz_raw(i,:) = cos(2*pi*i/T);
    test_Bz_raw(i,:) = sin(2*pi*i/T); test_Bt_raw(i,:) = cos(2*pi*i/T);
    %test_Bt_raw(i,:) = 1; test_Bz_raw(i,:) = 0;
    %test_Bt_raw(i,:) = 0; test_Bz_raw(i,:) = 1;
    %test_Bt_raw(i,:) = 2*sin(2*pi*i/T); test_Bz_raw(i,:) = sin(2*pi*i/T);

end



% signal format: n_time x n_probe (example: 1000 x 8 for real experiment)
% can subsitude test_signal with real signal in real experiment

%initialization
    alpha = zeros(length(t), probe_num);
    theta = zeros(length(t), probe_num);
    one_list = ones(1, probe_num);
    B = zeros(2,length(t),probe_num);
    

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

    %{
    initialization
    B = zeros(1,probe_num,2);
    mat = zeros(2,2,probe_num);
    
    mat = [cos(theta(i,:)) sin(theta(i,:)); -sin(theta(i,:)) cos(theta(i,:))];
    B = mat*[test_Bt_raw(i,:); test_Bz_raw(i,:)];
    Bt(i,:) = B(1,i,:);
    Bz(i,:) = B(2,i,:);
    %}
end

    f = figure;
   % for k = 1:probe_num
        k = 8; hold on;
        plot(t,test_Bt_raw(:,k),'o','color','r','markersize',3);
        plot(t,test_Bz_raw(:,k),'o','color','b','markersize',3);
        %legend('Bt_ raw','Bz_ raw');
        plot(t, Bt(:,k),'*','color','r','markersize',3);  
        plot(t, Bz(:,k),'*','color','b','markersize',3);
        legend('Bt','Bz');
        
   % end
    
    title('calibrated signal');
    hold off


