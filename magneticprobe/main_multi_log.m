%clear all
clearvars
folder_path=getenv('ts3u_path');%ts-3uまでのパス
%%%%%ここが各PCのパス
%環境変数を設定していない場合はパスを''内に全て記入する（使用しないパスは空白''で良い）
%pathname.ts3u='ts3u_path';%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier='fourier_path';%fourierのmd0（データックのショットが入ってる）までのpath
pathname.save='/Users/mgar/CAO_save'; %保存先

%%%%(1)spread sheetから ログのテーブルを取得してTに格納
%Github/test-open/getTS6log.mを使用
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);

IDXlist=[2911:2913 2925 2926 2927 2931 2933 2947:2950 2942 2943 2946];
IDX=IDXlist(1,15); %42
    date=T.date(IDX);
    shot=T.shot(IDX);
    TF_shot=T.TFoffset(IDX);
    offset_TF=isfinite(TF_shot);
    if isnan(T.EF_A_(IDX))%%NaNでないことを確認（ログが空白だとNaNになる）
        i_EF=150;
    else  %NaNなら150をとりあえず代入、記入されているときはその値を使う
        i_EF=T.EF_A_(IDX);
    end
% date=211223;
% TF_shot = 21 ;
% offset_TF = true;
% shot = 19; %28 37 38 43 44
% i_EF=150;

% ********* ALWAYS RUN THIS FUNCTION FIRST ***********
% parameters:(date,TF_shot,shot,offset_TF,offset_EF)
%[B_z,r_probe,z_probe,ch_dist,data] = get_B_z(200130,4,12,true,150);
%[B_z,r_probe,z_probe,ch_dist,data] = get_B_z(200202,7,6,true,0);
%[B_z,r_probe,z_probe,ch_dist,data,data_raw] = get_B_z(200718,27,30,true,70);
%[B_z,r_probe,z_probe,ch_dist,data,data_raw,shot_num] = get_B_z(date,TF_shot,shot,true,true);
[B_z,r_probe,z_probe,ch_dist,data,data_raw,shot_num] =get_B_z(date,TF_shot,shot,offset_TF,i_EF,folder_path);
%########## Raed oscilloscope (DL716) file ##########
% parameters:(date,shot,TF_shot,offset_TF)
%[low_n_data] = low_n_mode(date,shot,TF_shot,offset_TF);

% run this as well
% (ignore a column of broken probes)

B_z = B_z([2,3,4,6,7,8],2:end,:);
data = data([2,3,4,6,7,8],2:end,:);
z_probe = z_probe(2:end);
ch_dist = ch_dist([2,3,4,6,7,8],2:end);
r_probe = r_probe([2,3,4,6,7,8]);
%}

%## old ##
% B_z = B_z(:,2:end,:);
% data = data(:,2:end,:);
% z_probe = z_probe(2:end);
% ch_dist = ch_dist(:,2:end);

% ***** toroidal mode *****
%toroidal_mode_offset_new(low_n_data,shot_num,offset,smoothing,movemean,standarization);
%[low_n_signal] = toroidal_mode_offset_new(low_n_data,shot,true,false,true,false);

% ### plot contour figure ###
%parameters(low_n_signal,t_start,t_end)
%contour_low_n(low_n_signal,460,550);


% ************* PLOTTING FUNCTIONS *******************
% parameters:(B_z,ch_dist,start_time,end_time)
% plot_B_z_in_time(B_z,ch_dist,350,600);

% plot_psi_multi(B_z,r_probe,z_probe,461:1:475,true,true,true,false);
% plot_psi_multi(B_z,r_probe,z_probe,476:1:491,true,true,true,false);
% plot_psi_multi(B_z,r_probe,z_probe,492:1:507,true,true,true,false);
%plot_psi_multi(B_z,r_probe,z_probe,466:1:485,true,true,true,false);
%plot_psi_multi(B_z,r_probe,z_probe,461:1:480,true,true,true,false,shot);

plot_psi_xpoint(B_z,r_probe,z_probe,471:1:490,true,true,true,false,true,IDX);
plot_psi_xpoint(B_z,r_probe,z_probe,491:1:510,true,true,true,false,true,IDX);
%plot_psi_xpoint(B_z,r_probe,z_probe,501:1:520,true,true,true,false,true,IDX);

% parameters:(B_z,r_probe,z_probe,t,fitting,fill,fixed_Clayer,show_probe)
% plot_psi_at_t(B_z,r_probe,z_probe,476,true,true,true,true);

%parameters:(B_z,r_probe)
% plot_fitrate(B_z,r_probe,shot);



