%%%%%%%%%%%%%%%%%%%%%%%%
%200ch用新規pcbプローブのみでの磁気面（Bz）
%dtacqのshot番号を直接指定する場合
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%ここが各PCのパス
%【※コードを使用する前に】環境変数を設定しておくか、matlab内のコマンドからsetenv('パス名','アドレス')で指定してから動かす
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.save=getenv('savedata_path');%outputデータ保存先

pathname.rawdata38=getenv('rawdata038_path');%dtacq a038のrawdataの保管場所
pathname.woTFdata=getenv('woTFdata_path');%rawdata（TFoffset引いた）の保管場所

% setenv("rsGdrive","/Users/rsomeya/Library/CloudStorage/GoogleDrive-rsomeya2016@g.ecc.u-tokyo.ac.jp/マイドライブ")
pathname.rawdata=[getenv('rsGdrive') '/pcb'];%dtacqのrawdataの保管場所

% %直接入力の場合
dtacqlist=39;
shotlist=992;%991:992;%【input】dtacqの保存番号
tfshotlist=988;%0;
date = 230313;%【input】計測日
time = 470;%【input】表示時刻(us)
n_data=numel(shotlist);%計測データ数
EFlist = 150;%150;%【input】EF電流
TFlist = 0;
trange=430:590;%【input】計算時間範囲(ほぼ固定)
mesh_rz=50; %【input】rz方向のメッシュ数(ほぼ固定)
cmap = false;%【input】カラーマップの有無(ほぼ固定)

for i=1:n_data
    dtacq_num=dtacqlist;
    shot=shotlist(i);
    tfshot=tfshotlist;
    i_EF=EFlist;
    TF=TFlist;
    plot_psi200ch_at_t(time, date, dtacq_num, shot, tfshot, pathname,mesh_rz,i_EF,trange,TF,cmap);
end
