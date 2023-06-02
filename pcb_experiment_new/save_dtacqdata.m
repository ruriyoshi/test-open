%デジタイザ別の保存(2022/11/17)
%個別に環境変数a038_path, a039_path, a040_pathを設定する必要あり
clear all
%setenv('MDSPLUS_DIR','/usr/local/mdsplus');
%addpath(fullfile(getenv('MDSPLUS_DIR'), 'matlab'));
pathname.rawdata=getenv('rawdata_path');%保存先

dtacqlist=38;
shotlist=[3592:3594];%【input】dtacqの保存番号
%tfshotlist=zeros(size(shotlist));
tfshotlist = 3581*ones(size(shotlist));
date = 211223;%【input】計測日
n=numel(shotlist);%計測データ数

%RC係数読み込み

for i=1:n
    dtacq_num=dtacqlist;
    shot=shotlist(i);
    tfshot=tfshotlist(i);
    [rawdata]=getMDSdata(dtacq_num,shot,tfshot);%測定した生信号
    save(strcat(pathname.rawdata,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat'),'rawdata');
    if tfshot>0
        [rawdata0]=getMDSdata(dtacq_num,shot,0);
        save(strcat(pathname.rawdata,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot0.mat'),'rawdata0');
    end
end
