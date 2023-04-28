%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%グループ化した散布図をプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%各PCのパスを定義
run define_path.m

%------【input】---------------------------------------------------
date = 230315;%【input】実験日
begin_cal = 4;%【input】磁気面&フロー計算始めshot番号(実験ログD列)
end_cal = 35;%【input】磁気面&フロー計算終わりshot番号(実験ログD列)(0にするとbegin_cal以降の同日の全shot計算)
min_r = 12.5;%【input】ドップラープローブ計測点最小r座標[mm]
int_r = 2.5;%【input】ドップラープローブ計測点r方向間隔[mm]
min_z = -2.1;%【input】ドップラープローブ計測点最小z座標[mm](-2.1,2.1)
int_z = 4.2;%【input】ドップラープローブ計測点z方向間隔[mm](4.2)
ICCD.line = 'Ar';%【input】ドップラー発光ライン('Ar')
n_CH = 28;%【input】ドップラープローブファイバーCH数(28)
n_z = 1;%【input】ドップラープローブz方向データ数(数値)(1)
dtacq.num = 39;%【input】磁気プローブdtacq番号(39)
trange = 430:590;%【input】磁気プローブ計算時間範囲(430:590)
t_range = [477 483];%【input】計測時刻範囲(us)
r_range = [17.5 25];%【input】計測点範囲(cm)
ng_shotlist = [6 10 26 33 34];%【input】磁場失敗ショット番号
x_type = '|Br|';%【input】x軸の変数('Vz','Vr','|Vz|','|Vr|','Bz','Br','z','r','time')
y1_type = '|Vr|';%【input】y1軸の変数('Vz','Vr','|Vz|','|Vr|','Bz','Br','z','r','time')
y2_type = '|Vz|';%【input】y2軸の変数('Vz','Vr','|Vz|','|Vr|','Bz','Br','z','r','time',false)
g_type = 'time';%【input】グループ化基準('time','r')

%ドップラープローブ計測点配列を生成
mpoints = make_mpoints(n_CH,min_r,int_r,n_z,min_z,int_z);

%実験ログ読み取り
[exp_log,begin_row,end_row] = load_log(date);
if isempty(begin_row)
    return
end

figure('Position',[600 150 600 600])
if start_i <= end_row
    start_i = begin_row + begin_cal - 1;
    if end_cal == 0
        end_i = end_row;%begin_cal以降全部計算
    elseif end_cal < begin_cal
        error('end_cal must <= begin_cal.')
    elseif begin_row + end_cal - 1 <= end_row
        end_i = begin_row + end_cal - 1;%begin_calからend_calまで計算
    else
        error('end_cal must <= %d.', exp_log(end_row,4))
    end

    %配列を作成
    n_scatter = (end_i - start_i) * mpoints.n_r;%散布図プロット数
    sort = zeros(n_scatter,1);%範囲限定用
    V_z = zeros(n_scatter,1);%Vz
    V_r = zeros(n_scatter,1);%Vr
    T = zeros(n_scatter,1);%T
    B_z = zeros(n_scatter,1);%Bz
    B_r = zeros(n_scatter,1);%Br
    tp = zeros(n_scatter,1);%計測時刻
    zp = zeros(n_scatter,1);%計測点z座標
    rp = zeros(n_scatter,1);%計測点r座標
    for i = start_i:end_i
        step = i - start_i + 1;
        ICCD.shot = exp_log(i,4);%ショット番号
        a039shot = exp_log(i,8);%a039ショット番号
        a039tfshot = exp_log(i,9);%a039TFショット番号
        expval.PF1 = exp_log(i,11);%PF1電圧(kV)
        expval.PF2 = exp_log(i,14);%PF2電圧(kV)
        expval.TF = exp_log(i,18);%PF2電圧(kV)
        expval.EF = exp_log(i,23);%EF電流
        ICCD.trg = exp_log(i,42);%ICCDトリガ時間
        ICCD.exp_w = exp_log(i,43);%ICCD露光時間
        ICCD.gain = exp_log(i,44);%Andor gain
        time = round(ICCD.trg+ICCD.exp_w/2);%磁気面プロット時間
        if dtacq.num == 39
            dtacq.shot = a039shot;
            dtacq.tfshot = a039tfshot;
        end
        if time >= t_range(1) && time <= t_range(2) && isempty(intersect(ICCD.shot,ng_shotlist))
            %保存済みイオン温度、フローを読みこむ
            [V_i,absV,T_i,F,W,P,Lambda,Vx,Vy,ppoints,Angle] = load_ionvdist(date,ICCD,pathname);
            %磁場を読みこむ
            [grid2D,data2D,ok_z,ok_r] = load_pcb200ch(date,dtacq,pathname);
            for k = 1:mpoints.n_r
                if mpoints.r(k,1) >= r_range(1) && mpoints.r(k,1) <= r_range(2)
                    IDX_z = knnsearch(grid2D.zq(1,:).',mpoints.z(k,1)*1e-2);%grid2D.zqの中で最も計測z座標に近いセル番号を取得
                    IDX_r = knnsearch(grid2D.rq(:,1),mpoints.r(k,1)*1e-2);%grid2D.rqの中で最も計測r座標に近いセル番号を取得
                    IDX_t = (time - trange(1) + 1);
                    sort(mpoints.n_r*(step - 1) + k,1) = 1;
                    V_z(mpoints.n_r*(step - 1) + k,1) = V_i(k,1);%Vz
                    V_r(mpoints.n_r*(step - 1) + k,1) = V_i(k,2);%Vr
                    T(mpoints.n_r*(step - 1) + k,1) = T_i(k,1);%イオン温度
                    B_z(mpoints.n_r*(step - 1) + k,1) = data2D.Bz(IDX_r,IDX_z,IDX_t)*1e3;%Bz[mT]
                    B_r(mpoints.n_r*(step - 1) + k,1) = data2D.Br(IDX_r,IDX_z,IDX_t)*1e3;%Br[mT]
                    tp(mpoints.n_r*(step - 1) + k,1) = trange(IDX_t);%計測時刻[us]
                    zp(mpoints.n_r*(step - 1) + k,1) = mpoints.z(k,1);%z座標[cm]
                    rp(mpoints.n_r*(step - 1) + k,1) = mpoints.r(k,1);%r座標[cm]
                end
            end
            % else
            %     V_z(mpoints.n_r*(step - 1) + k,1) = NaN;
            %     V_r(mpoints.n_r*(step - 1) + k,1) = NaN;
            %     T(mpoints.n_r*(step - 1) + k,1) = NaN;
            %     B_z(mpoints.n_r*(step - 1) + k,1) = NaN;
            %     B_r(mpoints.n_r*(step - 1) + k,1) = NaN;
            %     tp(mpoints.n_r*(step - 1) + k,1) = NaN;
            %     zp(mpoints.n_r*(step - 1) + k,1) = NaN;
            %     rp(mpoints.n_r*(step - 1) + k,1) = NaN;
        end
    end
    indices = find(sort==0);
    V_z(indices) = [];
    V_r(indices) = [];
    T(indices) = [];
    B_z(indices) = [];
    B_r(indices) = [];
    tp(indices) = [];
    zp(indices) = [];
    rp(indices) = [];

    %グループ分けされた散布図を作成
    switch x_type
        case 'Vz'
            x_val = V_z;
            x_name = 'V_z [km/s]';
        case 'Vr'
            x_val = V_r;
            x_name = 'V_r [km/s]';
        case '|Vz|'
            x_val = abs(V_z);
            x_name = '|V_z| [km/s]';
        case '|Vr|'
            x_val = abs(V_r);
            x_name = '|V_r| [km/s]';
        case 'T'
            x_val = T;
            x_name = 'T [eV]';
        case 'Bz'
            x_val = B_z;
            x_name = 'B_z [mT]';
        case 'Br'
            x_val = B_r;
            x_name = 'B_r [mT]';
        case '|Bz|'
            x_val = abs(B_z);
            x_name = '|B_z| [mT]';
        case '|Br|'
            x_val = abs(B_r);
            x_name = '|B_r| [mT]';
        case 'time'
            x_val = tp;
            x_name = 'time [us]';
        case 'z'
            x_val = zp;
            x_name = 'z [cm]';
        case 'r'
            x_val = rp;
            x_name = 'r [cm]';
        otherwise
            warning('Input error in x_type.')%x_typeの入力エラー
            return;
    end
    switch y1_type
        case 'Vz'
            y1_val = V_z;
            y1_name = 'V_z [km/s]';
        case 'Vr'
            y1_val = V_r;
            y1_name = 'V_r [km/s]';
        case '|Vz|'
            y1_val = abs(V_z);
            y1_name = '|V_z| [km/s]';
        case '|Vr|'
            y1_val = abs(V_r);
            y1_name = '|V_r| [km/s]';
        case 'T'
            y1_val = T;
            y1_name = 'T [eV]';
        case 'Bz'
            y1_val = B_z;
            y1_name = 'B_z [mT]';
        case 'Br'
            y1_val = B_r;
            y1_name = 'B_r [mT]';
        case '|Bz|'
            y1_val = abs(B_z);
            y1_name = '|B_z| [mT]';
        case '|Br|'
            y1_val = abs(B_r);
            y1_name = '|B_r| [mT]';
        case 'time'
            y1_val = tp;
            y1_name = 'time [us]';
        case 'z'
            y1_val = zp;
            y1_name = 'z [cm]';
        case 'r'
            y1_val = rp;
            y1_name = 'r [cm]';
        otherwise
            warning('Input error in y1_type.')%y1_typeの入力エラー
            return;
    end
    switch g_type
        case 'time'
            g_val = tp;
            g_name = 'time [us]';
        case 'z'
            g_val = zp;
            g_name = 'z [cm]';
        case 'r'
            g_val = rp;
            g_name = 'r [cm]';
        otherwise
            warning('Input error in g_type.')%g_typeの入力エラー
            return;
    end
    n_group = length(unique(g_val, 'rows'));
    clr = hsv(n_group);
    gs1 = gscatter(x_val,y1_val,g_val,clr,'o','filled');
    hold on
    if y2_type
        switch y2_type
            case 'Vz'
                y2_val = V_z;
                y2_name = 'V_z [km/s]';
            case 'Vr'
                y2_val = V_r;
                y2_name = 'V_r [km/s]';
            case '|Vz|'
                y2_val = abs(V_z);
                y2_name = '|V_z| [km/s]';
            case '|Vr|'
                y2_val = abs(V_r);
                y2_name = '|V_r| [km/s]';
            case 'T'
                y2_val = T;
                y2_name = 'T [eV]';
            case 'Bz'
                y2_val = B_z;
                y2_name = 'B_z [mT]';
            case 'Br'
                y2_val = B_r;
                y2_name = 'B_r [mT]';
            case '|Bz|'
                y2_val = abs(B_z);
                y2_name = '|B_z| [mT]';
            case '|Br|'
                y2_val = abs(B_r);
                y2_name = '|B_r| [mT]';
            case 'time'
                y2_val = tp;
                y2_name = 'time [us]';
            case 'z'
                y2_val = zp;
                y2_name = 'z [cm]';
            case 'r'
                y2_val = rp;
                y2_name = 'r [cm]';
            otherwise
                warning('Input error in y2_type.')%y2_typeの入力エラー
                return;
        end
        gs2 = gscatter(x_val,y2_val,g_val,clr,'^');
        y_name = [y1_name,'  ',y2_name];
        title_name = ['(',y1_name,', ',y2_name,') VS ',x_name,' grouped by ',g_name];
    else
        y_name = y1_name;
        title_name = [y1_name,' VS ',x_name,' grouped by ',g_name];
    end
    title(title_name)
    % yline(0)%直線y = 0を挿入
    xlabel(x_name)
    ylabel(y_name)
    lgd = legend('Location','eastoutside');
    title(lgd,g_name)
    legend
    ax = gca;
    ax.FontSize = 14;
    hold off
else
    error('begin_cal must <= %d.', exp_log(end_row,4))
end
