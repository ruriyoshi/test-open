%%%%%%%%%%%pcb計測の二次元データdata2D（Bz,Br，Et，Jt）とグリッドgrid2D（zq,rq）が出てくる。
%%%% input （日付、データックのショット番号、差し引きショット、見たい時間、nはグリッドの個数、EFの電流）
%%%%%%↑　これはtestの中にある。psisave.mlx内のfunction [date, shot, TF_shot, offset_TF,i_EF,start,Doppler_t] = getinput(T,IDX)を関数化して呼び出すのが楽
function [grid2D, data2D] = pcbdata(date, d_tacq,d_tacqTF,trange, coeff, n,EF)
%[rawdata]=getvalue(d_tacq,d_tacqTF); % rawdata
pathname.rawdata='C:\Users\uswk0\OneDrive\デスクトップ\data\rawdata_a038\'; %rawdataの保管場所（各自変更）
filename=strcat(pathname.rawdata,'rawdata_dtacq',num2str(d_tacq),'.mat');
load(filename,'rawdata');

if numel(rawdata)< 500
    grid2D=NaN;
    data2D=NaN;
    return
end
if exist(coeff)==0
    load('rc_coeff2020.mat')
end     
[ok, bz, rpos, zpos,p_ch] = getpcbbz(rawdata, coeff,date);

% チャンネルごとの生信号のプロット
% bz_s=smoothdata(bz,1);
% r = 5;
% col1 = 12;
% col2 = 13;
% y_upper_lim = 0.1;
% y_lower_lim = -0.1;
% t_start=300;
% t_end=600;
% r_ch=col1+col2;
% plotbzsignal(y_upper_lim, col2, col1, t_end, p_ch, y_lower_lim, t_start, bz_s, ok, r_ch, r);
% clear r col1 col2 y_upper_lim y_lower_lim t t_start t_end

%zq*rqの座標を生成
[zq,rq]=meshgrid(linspace(min(zpos),max(zpos),n),linspace(min(rpos),max(rpos),n));
%構造体にして情報を入れ込む
grid2D=struct('zq',zq,'rq',rq);
clear zq rq
%data2D（Bz,Br，Et，Jt）を計算してる
data2D = data2Dcalc(EF, grid2D, n, trange, rpos, zpos, bz, ok);
end
