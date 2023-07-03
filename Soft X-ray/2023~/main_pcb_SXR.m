% 土井プローブ(Bz100ch)の磁気面に，SXR画像を重ねて表示します．
clear

pathname.NIFS=getenv('NIFS_path');% 192.168.1.111
pathname.rawdata=getenv('rawdata_path');% dtacqのrawdataの保管場所
addpath(genpath('/Users/yuleo/Documents/GitHub/test-open'));
% 実験オペレーションの取得
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
pat = 230315;
date = pat;
T=searchlog(T,node,pat);
% IDXlist = Indexlistは，表示させたいshot番号を意味します．これを用いて，その他の殆どのパラメータをスプレッドシートから自動で読みよることを試みます．
IDXlist = 33;% [33,35:37]
n_data=numel(IDXlist);% n_dataは計測データ数を意味します．
shotlist=T.a039(IDXlist);
tfshotlist=T.a039_TF(IDXlist);
EFlist=T.EF_A_(IDXlist);
TFlist=T.TF_kV_(IDXlist);
dtacqlist=39.*ones(n_data,1);
startlist = T.SXRStart(IDXlist);
intervallist = T.SXRInterval(IDXlist);
% 
% 直接入力の場合【注意】全て同じサイズの行列になるように記入
% dtacqlist=39;
% shotlist=1118;%【input】実験ログのa039の番号
% tfshotlist=1106;%【input】実験ログのa039_TFの番号
% date = 230315;%【input】計測日
% n_data=numel(shotlistかかかかあｋ);%計測データ数
% EFlist = 150;%【input】EF電流

trange=440:500;%【input】計算時間範囲
n=50; %【input】rz方向のメッシュ数

t = 475;
show_xpoint = false;
show_localmax = false;
% start = 450;
% interval = 5;
save = true;
FilterProcessFlag = false;
NonLinerFunctionFlag = false;

for i=1:n_data
    dtacq_num=dtacqlist(i);
    shot=shotlist(i);
    tfshot=tfshotlist(i);
    i_EF=EFlist(i);
    TF=TFlist(i);
    start = startlist(i);
    interval = intervallist(i);
%     TF=4;
%     plot_psi200ch(date, dtacq_num, shot, tfshot, pathname,n,i_EF,trange,TF); 
    [grid2D,data2D] = process_PCBdata(date, dtacq_num, shot, tfshot, pathname, n,i_EF,trange);
    shot_SXR = IDXlist(i);
%     shot_SXR = 14;
    SXRImageFilePath = strcat('G:/My Drive/X-ray/Data/SXROUT/',num2str(date),'/shot',num2str(shot_SXR,'%03i'),'.tif');
    [EE_high,EE_low] = plot_SXR_at_t(grid2D,data2D,date,shot_SXR,t,show_xpoint,show_localmax,start,interval,save,SXRImageFilePath,FilterProcessFlag,NonLinerFunctionFlag);
%    plot_SXR_multi(grid2D,data2D,date,shot_SXR,show_xpoint,show_localmax,start,interval,save,SXRfilename,filter,NL);
%     Brec = clc_Breconnection(grid2D,data2D);
end
% figure;plot(data2D.trange,Brec);