%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ショット番号、撮影パラメータなどを実験ログから自動取得して
%ドップラープローブによるイオン速度分布関数、温度、フローとその瞬間の磁気面をプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

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

%------【input】---------------------------------------------------
date = 230310;%【input】実験日
begin_cal = 1;%【input】磁気面&フロー計算始めshot番号(実験ログD列)
end_cal = 0;%【input】磁気面&フロー計算終わりshot番号(実験ログD列)(0にするとbegin_cal以降の同日の全shot計算)
min_r = 12.5;%【input】ドップラープローブ計測点最小r座標[mm]
int_r = 2.5;%【input】ドップラープローブ計測点r方向間隔[mm]
min_z = 2.1;%【input】ドップラープローブ計測点最小z座標[mm](-2.1,2.1)
int_z = 4.2;%【input】ドップラープローブ計測点z方向間隔[mm](4.2)
ICCD.line = 'Ar';%【input】ドップラー発光ライン('Ar')
n_CH = 28;%【input】ドップラープローブファイバーCH数(28)
n_z = 1;%【input】ドップラープローブz方向データ数(数値)(1)
%-----------------------詳細設定【input】----------------------------
cal_vdist = true;%【input】速度分布を計算(true,false)
plot_spectra = false;%【input】スペクトルをプロット(true,false)
plot_analisis = false;%【input】逆変換解析をプロット(true,false)
plot_vdist = false;%【input】速度分布をプロット(true,false)
plot_compare = false;%【input】再構成比較をプロット(true,false)
plot_flow = true;%【input】流速をプロット(true,false)
plot_psi = false;%【input】磁気面をプロット(true,false)
overlay_plot = false;%【input】流速と磁気面を重ねる(true,false)

save_fit = false;%【input】ガウスフィッティングpngを保存(true,false)
save_fig = true;%【input】速度分布、フローpngを保存(true,false)

save_vdist = true;%【input】速度分布データを保存(true,false)
load_vdist = false;%【input】速度分布データを読み込む(true,false)

plot_type = 'contour';%【input】速度分布プロット種類('contour','surf')
Ti_type = 'dispersion';%【input】イオン温度計算法('dispersion')

show_offset = false;%【input】分光offsetを表示(true,false)
inversion_method = 5;%【input】速度分布逆変換手法(1~6)
factor = 0.1;%【input】イオンフロー矢印サイズ(数値:0.1など)
dtacq.num = 39;%【input】磁気プローブdtacq番号(39)
mesh_rz = 50;%【input】磁気プローブrz方向のメッシュ数(50)
trange = 430:590;%【input】磁気プローブ計算時間範囲(430:590)
%------------------------------------------------------------------

%ドップラープローブ計測点配列を生成
mpoints = make_mpoints(n_CH,min_r,int_r,n_z,min_z,int_z);

%実験ログ読み取り
[exp_log,begin_row,end_row] = load_log(date);
if isempty(begin_row)
    return
end

%--------磁気面&フローを計算------
start_i = begin_row + begin_cal - 1;
if start_i <= end_row
    if end_cal == 0
        end_i = end_row;%begin_cal以降全部計算
    elseif end_cal < begin_cal
        error('end_cal must <= begin_cal.')
    elseif begin_row + end_cal - 1 <= end_row
        end_i = begin_row + end_cal - 1;%begin_calからend_calまで計算
    else
        error('end_cal must <= %d.', exp_log(end_row,4))
    end
    for i = start_i:end_i
        ICCD.shot = exp_log(i,4);%ショット番号
        a039shot = exp_log(i,8);%a039ショット番号
        a039tfshot = exp_log(i,9);%a039TFショット番号
        expval.PF1 = exp_log(i,11);%PF1電圧(kv)
        expval.PF2 = exp_log(i,14);%PF2電圧(kv)
        expval.TF = exp_log(i,18);%PF2電圧(kv)
        expval.EF = exp_log(i,23);%EF電流
        ICCD.trg = exp_log(i,42);%ICCDトリガ時間
        ICCD.exp_w = exp_log(i,43);%ICCD露光時間
        ICCD.gain = exp_log(i,44);%Andor gain
        time = round(ICCD.trg+ICCD.exp_w/2);%磁気面プロット時間
        if dtacq.num == 39
            dtacq.shot = a039shot;
            dtacq.tfshot = a039tfshot;
        end
        if cal_vdist
            %イオン速度分布を計算
            [V_i,absV,T_i] = cal_ionvdist(date,expval,ICCD,mpoints,pathname,show_offset,plot_spectra, ...
                inversion_method,plot_analisis,plot_vdist,plot_type,save_fig,plot_compare,save_vdist,Ti_type);
        elseif load_vdist
            %保存済みイオン温度、フローを読み取り
            [V_i,absV,T_i,F,W,P,Lambda,Vx,Vy,ppoints,Angle] = load_ionvdist(date,ICCD,pathname);
            if plot_vdist
                plot_ionvdist(Vx,Vy,F,date,expval,ICCD,pathname,mpoints,ppoints,plot_type,save_fig)
            end
            if plot_compare
                plot_inversion_compare(F,W,P,Lambda,mpoints,Angle,ICCD)
            end
        else
            V_i = char.empty;
            absV = char.empty;
            T_i = char.empty;
        end
        %磁気面をプロット
        if plot_psi
            plot_psi200ch_at_t(time,date,dtacq,pathname,mesh_rz,expval,trange,false);
        end
        %イオン温度、フローをプロット
        if plot_flow
            if not(isempty(V_i))
                if plot_psi
                    plot_ionflow(V_i,absV,T_i,date,expval,ICCD,pathname,factor,mpoints,overlay_plot,save_fig,'ionvdist')
                else
                    plot_ionflow(V_i,absV,T_i,date,expval,ICCD,pathname,factor,mpoints,false,save_fig,'ionvdist')
                end
            end
        end
    end
else
    error('begin_cal must <= %d.', exp_log(end_row,4))
end
