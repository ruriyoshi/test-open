function [ch1,ch2,ch3,ch4,ch5,ch6,ch7,ch8] = get_oscillo(date,TF_only,shot,offset_TF,offset_EF)
%GET_IP 
% input:
%   string: directory of the interested raw data rgw file
%   string: directory of TF discharge only raw data rgw file
% output:
%   3d array of double: B_z (r,z,t), offsetted at zero and smoothed

% read file

if shot < 10
    data_dir = ['/Users/keisukemiki/koala_mnt/experiment/results/ts-3u/',num2str(date),'/',num2str(date),'00',num2str(shot),'.rgw'];
elseif shot < 100
    data_dir = ['/Users/keisukemiki/koala_mnt/experiment/results/ts-3u/',num2str(date),'/',num2str(date),'0',num2str(shot),'.rgw'];
elseif shot < 1000
    data_dir = ['/Users/keisukemiki/koala_mnt/experiment/results/ts-3u/',num2str(date),'/',num2str(date),num2str(shot),'.rgw'];
else
    disp('More than 999 shots! You need some rest!!!')
end
data=dlmread(data_dir,'\t',1,1);

% set data
ch1=data(:,2);
ch2=data(:,3);
ch3=data(:,4);
ch4=data(:,5);
ch5=data(:,6);
ch6=data(:,7);
ch7=data(:,8);
ch8=data(:,9);

%ip_current = (ch1*0.1);
%cep_test = (ch4*-782.5873*0.8222);
%cep1 = (ch5);%*7.1644e+02);
%cep2 = (ch6);%*2.3020e+01);
%pf1 = (ch7);%*-1.5024e+03);
%pf2 = (ch8);%*-2.3876e+03);
end