%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ショット番号、撮影パラメータを手入力して
%ドップラープローブによるイオン温度、フローをプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%ここが各PCのパス
%【※コードを使用する前に】環境変数を設定しておくか、matlab内のコマンドからsetenv('パス名','アドレス')で指定してから動かす
setenv("NIFS_path","/Volumes/experiment/results")
setenv("rsGdrive","/Users/rsomeya/Library/CloudStorage/GoogleDrive-rsomeya2016@g.ecc.u-tokyo.ac.jp/マイドライブ/lab")
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.save=[getenv('rsGdrive') '/save'];%outputデータ保存先
pathname.rawdata38=getenv('rawdata038_path');%dtacq a038のrawdataの保管場所
pathname.woTFdata=getenv('woTFdata_path');%rawdata（TFoffset引いた）の保管場所
pathname.fig=[getenv('rsGdrive') '/figure'];%figure保存先
pathname.mat=[getenv('rsGdrive') '/mat'];%figure保存先
pathname.rawdata=[pathname.mat,'/pcb'];%dtacqのrawdataの保管場所
pathname.flowdata=[pathname.mat,'/ionflow'];%流速データの保管場所
pathname.vdistdata=[pathname.mat,'/ionvdist'];%速度分布データの保管場所

%------【input】-------
date = 230315;%【input】実験日
expval.PF1 = 39;%PF1電圧(kv)
expval.PF2 = 39;%PF2電圧(kv)
expval.TF = 4;%PF2電圧(kv)
expval.EF = 150;%EF電流
ICCD.shot = 5;%【input】ショット番号
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
%------詳細設定【input】-------
cal_flow = true;%【input】流速を計算(true,false)
show_offset = false;%【input】分光offsetを表示(true,false)
plot_fit = true;%【input】ガウスフィッティングを表示(true,false)
save_fit = false;%【input】ガウスフィッティングpngを保存(true,false)
save_flow = false;%【input】流速データを保存(true,false)
load_flow = false;%【input】流速データを読み込む(true,false)
plot_flow = false;%【input】流速をプロット(true,false)
save_fig = false;%【input】流速pngを保存(true,false)
factor = 0.1;%【input】イオンフロー矢印サイズ(数値:0.1など)

%計測点配列を生成
mpoints = make_mpoints(n_CH,min_r,int_r,n_z,min_z,int_z);

%イオン温度、フローを計算、プロット
if cal_flow
    %イオン温度、フローを計算
    [V_i,absV,T_i] = cal_ionflow(date,ICCD,mpoints,pathname,show_offset,plot_fit,save_fit,save_flow);
elseif load_flow
    %保存済みイオン温度、フローを読み取り
    [V_i,absV,T_i] = load_ionflow(date,ICCD,pathname);
end

if plot_flow
    if isempty(V_i)
        return
    else
        plot_ionflow(V_i,absV,T_i,date,expval,ICCD,pathname,factor,mpoints,false,save_fig,'ionflow')
    end
end