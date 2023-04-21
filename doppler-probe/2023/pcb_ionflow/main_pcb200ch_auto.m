%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dtacqショット番号を実験ログから
%自動取得して磁気面をプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
date = 230313;%【input】実験日(yymmdd)
begin_cal = 5;%【input】計算始めshot番号(実験ログD列)
end_cal = 5;%【input】計算終わりshot番号(実験ログD列)(0にするとbegin_cal以降同日の全shot計算)
Nplot = 16;%【input】磁気面subplot枚数(16以下の4の倍数or3以下)
t_start = 460;%【input】磁気面subplot開始時間(us)
dt = 2;%【input】磁気面subplot時間間隔(us)
cmap = false;%【input】磁気面カラーマップ('psi','Bz','Bt','Jt','Et',false)
cbar = true;%【input】カラーバー(true,false)
dtacq_num = 39;%【input】磁気プローブdtacq番号(ほぼ固定)
mesh_rz = 50;%【input】磁気プローブrz方向のメッシュ数(ほぼ固定)
trange = 430:590;%【input】磁気プローブ計算時間範囲(ほぼ固定)

%実験ログ読み取り
[exp_log,begin_row,end_row] = load_log(date);
if isempty(begin_row)
    return
end

%--------磁気面を計算------
start_i = begin_row + begin_cal - 1;
if start_i <= end_row
    if end_cal == 0
        end_i = end_row;%begin_cal以降全部計算
    elseif end_cal < begin_cal
        warning('end_cal must <= begin_cal.')
        return
    elseif begin_row + end_cal - 1 <= end_row
        end_i = begin_row + end_cal - 1;
    else
        warning('end_cal must <= %d.', exp_log(end_row,4))
        return
    end
    for i = start_i:end_i
        % shot = exp_log(i,4);%ショット番号
        a039shot = exp_log(i,8);%a039ショット番号
        a039tfshot = exp_log(i,9);%a039TFショット番号
        i_EF = exp_log(i,23);%EF電流
        % ICCD.trg = exp_log(i,42);%ICCDトリガ時間
        % ICCD.exp_w = exp_log(i,43);%ICCD露光時間
        % ICCD.gain = exp_log(i,44);%Andor gain
        if dtacq_num == 39
            dtacq_shot = a039shot;
            dtacq_tfshot = a039tfshot;
        end
        plot_psi200ch_multi(Nplot,t_start,dt,date,dtacq_num,dtacq_shot,dtacq_tfshot,pathname,mesh_rz,i_EF,trange,cmap,cbar);
    end
else
    warning('begin_cal must <= %d.', exp_log(end_row,4))
    return
end
