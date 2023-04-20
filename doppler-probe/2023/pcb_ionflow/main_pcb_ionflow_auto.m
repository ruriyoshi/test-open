%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ショット番号、撮影パラメータなどを実験ログから自動取得して
%ドップラープローブによるイオン温度、フローとその瞬間の磁気面をプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = main_pcb_ionflow_auto(show_offset,plot_fit,plot_flow,save_flow,save_fig,plot_psi)
%offsetを表示/ガウスフィッティングを表示/流速をプロット/流速データを保存/流速figを保存/磁気面をプロット
%全てtrue or false

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
date = 230313;%【input】実験日
begin_cal = 3;%【input】磁気面&フロー計算始めshot番号(実験ログD列)
end_cal = 3;%【input】磁気面&フロー計算終わりshot番号(実験ログD列)(0にするとbegin_cal以降の同日の全shot計算)
min_r = 12.5;%【input】ドップラープローブ計測点最小r座標[mm]
int_r = 2.5;%【input】ドップラープローブ計測点r方向間隔[mm]
min_z = 2.1;%【input】ドップラープローブ計測点最小z座標[mm]
int_z = 4.2;%【input】ドップラープローブ計測点z方向間隔[mm]
gas = 'Ar';%【input】ガス種('Ar')
NofCH = 28;%【input】ファイバーCH数(28)
nz = 1;%【input】z方向データ数(数値)(1)
factor = 0.05;%【input】矢印サイズ(数値:0.05など)
dtacq_num = 39;%【input】磁気プローブdtacq番号
mesh_rz = 50;%【input】磁気プローブrz方向のメッシュ数(ほぼ固定)
trange = 430:590;%【input】磁気プローブ計算時間範囲(ほぼ固定)

%計測点配列を生成
[r_measured,z_measured] = make_mpoints(NofCH,min_r,int_r,nz,min_z,int_z);

%--------------実験ログ読み取り---------------
%時短用(固定値で良い)
load_s = 5300;%実験ログ読み始め行番号(230309~)
load_f = 10000;%実験ログ読み終わり行番号(固定)
%実験ログがcdになければ最新の実験ログをwebから取得
if exist('exp_log.xlsx','file')
    FileInfo = dir('exp_log.xlsx');
    DateNum = FileInfo.datenum;
    FileDate = datetime(DateNum, 'ConvertFrom', 'datenum', 'Format', 'yyMMdd');
    if str2double(string(FileDate)) > date
        disp(append(string(FileDate), '更新の実験ログを使用します。'))
    else
        run save_log.m
        disp('最新の実験ログ(exp_log.xlsx)を保存しました。')
    end
else
    run save_log.m
    disp('最新の実験ログ(exp_log.xlsx)を保存しました。')
end
%実験ログ中の実験日に対応する範囲を特定
exp_log = readmatrix('exp_log.xlsx','Sheet','log','Range', ['A' num2str(load_s) ':AR' num2str(load_f)]);
[n_row,n_col] = size(exp_log);
begin_row = find(exp_log(:,3) == date);%実験日の最初のshotの行番号を取得
if isempty(begin_row)
    warning('実験日が実験ログ中に存在しません。')
    return
end
end_row = begin_row;
while end_row<n_row && isnan(exp_log(end_row+1,3)) && exp_log(end_row+1,4)%日付がNaN&&shot番号が記入済=実験日のshot
    end_row = end_row+1;
end%実験日の最後のshotの行番号を取得

%--------磁気面&フローを計算------
start_i = begin_row + begin_cal - 1;
if start_i <= end_row
    if end_cal == 0
        end_i = end_row;%begin_cal以降全部計算
    elseif end_cal < begin_cal
        warning('end_cal must <= begin_cal.')
        return
    elseif begin_row + end_cal - 1 <= end_row
        end_i = begin_row + end_cal - 1;%begin_calからend_calまで計算
    else
        warning('end_cal must <= %d.', exp_log(end_row,4))
        return
    end
    for i = start_i:end_i
        shot = exp_log(i,4);%ショット番号
        a039shot = exp_log(i,8);%a039ショット番号
        a039tfshot = exp_log(i,9);%a039TFショット番号
        i_EF = exp_log(i,23);%EF電流
        trg = exp_log(i,42);%ICCDトリガ時間
        exp_w = exp_log(i,43);%ICCD露光時間
        gain = exp_log(i,44);%Andor gain
        if dtacq_num == 39
            dtacq_shot = a039shot;
            dtacq_tfshot = a039tfshot;
        end
        if plot_psi
            plot_psi200ch_at_t(round(trg+exp_w/2), date, dtacq_num, dtacq_shot, dtacq_tfshot, pathname,mesh_rz,i_EF,trange,false);
        end
        plot_ionflow(date,shot,trg,exp_w,gain,gas,NofCH,nz,show_offset,plot_fit,plot_flow,save_flow,save_fig,factor,r_measured,z_measured);
    end
else
    warning('begin_cal must <= %d.', exp_log(end_row,4))
    return
end

