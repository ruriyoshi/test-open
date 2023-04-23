%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ショット番号、撮影パラメータを手入力して
%ドップラープローブによるイオン温度、フローをプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = main_ionflow_manual(cal_flow,save_flow,load_flow,plot_flow,save_fig)
%流速を計算/流速データを保存/流速データを読み込み/流速をプロット/流速figを保存
%全てtrue or false
%実行例)main_ionflow_manual(true,true,false,true,false)流速データを保存してプロット
%実行例)main_ionflow_manual(false,false,true,true,false)保存済み流速データからプロット

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
pathname.flowdata=[getenv('rsGdrive') '/ionflow'];%流速データの保管場所

%------【input】-------
date = 230313;%【input】実験日
ICCD.shot = 4;%【input】ショット番号
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
factor = 0.05;%【input】イオンフロー矢印サイズ(数値:0.05など)
show_offset = false;%【input】分光offsetを表示(true,false)
plot_fit = false;%【input】ガウスフィッティングを表示(true,false)

%計測点配列を生成
mpoints = make_mpoints(n_CH,min_r,int_r,n_z,min_z,int_z);

%イオン温度、フローを計算、プロット
if cal_flow
    %イオン温度、フローを計算
    [V_i,absV,T_i] = cal_ionflow(date,ICCD,mpoints,pathname,show_offset,plot_fit,save_flow);
elseif load_flow
    %保存済みイオン温度、フローを読み取り
    [V_i,absV,T_i] = load_ionflow(date,ICCD,pathname);
end

if plot_flow
    if isempty(V_i)
        return
    else
        plot_ionflow(V_i,absV,T_i,date,ICCD,factor,mpoints,false,save_fig)
    end
end