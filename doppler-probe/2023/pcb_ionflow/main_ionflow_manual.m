%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ショット番号、撮影パラメータを手入力して
%ドップラープローブによるイオン温度、フローをプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = main_ionflow_manual(show_offset,plot_fit,save_flow,plot_flow,save_fig)
%分光offsetを表示/ガウスフィッティングを表示/流速データを保存/流速をプロット/流速figを保存
%全てtrue or false
%実行例)main_ionflow_manual(false,false,false,true,false)

%%%%%ここが各PCのパス
%【※コードを使用する前に】環境変数を設定しておくか、matlab内のコマンドからsetenv('パス名','アドレス')で指定してから動かす
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
% setenv("NIFS_path","/Volumes/experiment/results")
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.save=getenv('savedata_path');%outputデータ保存先
pathname.rawdata38=getenv('rawdata038_path');%dtacq a038のrawdataの保管場所
pathname.woTFdata=getenv('woTFdata_path');%rawdata（TFoffset引いた）の保管場所
% setenv("rsGdrive","/Users/rsomeya/Library/CloudStorage/GoogleDrive-rsomeya2016@g.ecc.u-tokyo.ac.jp/マイドライブ")
pathname.rawdata=[getenv('rsGdrive') '/pcb'];%dtacqのrawdataの保管場所

%------【input】-------
date = 230313;%【input】実験日
ICCD.shot = 3;%【input】ショット番号
ICCD.trg = 474;%【input】ICCDトリガ時間
ICCD.exp_w = 2;%【input】ICCD露光時間
ICCD.gain = 4095;%【input】ICCD gain
ICCD.line = 'Ar';%【input】ドップラー発光ライン('Ar')
min_r = 12.5;%【input】ドップラープローブ計測点最小r座標[cm]
int_r = 2.5;%【input】ドップラープローブ計測点r方向間隔[cm]
min_z = 2.1;%【input】ドップラープローブ計測点最小z座標[cm]
int_z = 4.2;%【input】ドップラープローブ計測点z方向間隔[cm]
n_CH = 28;%【input】ドップラープローブファイバーCH数(28)
n_z = 1;%【input】ドップラープローブz方向データ数(数値)(1)
factor = 0.05;%【input】イオンフロー矢印サイズ(数値:0.05など)

%計測点配列を生成
mpoints = make_mpoints(n_CH,min_r,int_r,n_z,min_z,int_z);

%イオン温度、フローを計算、プロット
[V_i,absV,T_i] = cal_ionflow(date,ICCD,mpoints,pathname,show_offset,plot_fit,save_flow);
if plot_flow
    plot_ionflow(V_i,absV,T_i,date,ICCD,factor,mpoints,false,save_fig)
end
end