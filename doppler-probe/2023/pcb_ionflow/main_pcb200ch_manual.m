%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dtacqショット番号を手入力して
%磁気面をプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%------【input】-------
date = 230313;%【input】実験日(yymmdd)
shotlist=991;%991:992;%【input】dtacqの保存番号
tfshotlist=988;%【input】dtacqのtfonly保存番号
EFlist = 150;%150;%【input】EF電流
Nplot = 16;%【input】磁気面subplot枚数(16以下の4の倍数or3以下)
t_start = 460;%【input】磁気面subplot開始時間(us)
dt = 2;%【input】磁気面subplot時間間隔(us)
cmap = false;%【input】磁気面カラーマップ('psi','Bz','Bt','Jt','Et',false)
cbar = true;%【input】カラーバー(true,false)
%------詳細設定【input】-------
dtacq_num = 39;%【input】磁気プローブdtacq番号(ほぼ固定)
mesh_rz = 50;%【input】磁気プローブrz方向のメッシュ数(ほぼ固定)
trange = 430:590;%【input】磁気プローブ計算時間範囲(ほぼ固定)

n_data=numel(shotlist);%計測データ数
for i=1:n_data
    dtacq_shot=shotlist(i);
    dtacq_tfshot=tfshotlist;
    i_EF=EFlist;
    plot_psi200ch_multi(Nplot,t_start,dt,date,dtacq_num,dtacq_shot,dtacq_tfshot,pathname,mesh_rz,i_EF,trange,cmap,cbar);
end
