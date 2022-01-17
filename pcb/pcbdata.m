function [grid2D, data2D] = pcbdata(date, shot,tfshot,trange, coeff, n,EF)
[rawdata]=getvalue(shot,tfshot); % rawdata
if numel(rawdata)< 500
    grid2D=NaN;
    data2D=NaN;
    return
end
if exist(coeff)==0
    load('rc_coeff2020.mat')
end     
[ok, bz, rpos, zpos,p_ch] = getpcbbz(rawdata, coeff,date);
% チャンネルごとの信号
%bz=smoothdata(bz,1);
% r = 5;
% col1 = 12;
% col2 = 13;
% y_upper_lim = 0.1;
% y_lower_lim = -0.1;
% t_start=300;
% t_end=600;
% r_ch=col1+col2;
% plotbzsignal(y_upper_lim, col2, col1, t_end, p_ch, y_lower_lim, t_start, bz, ok, r_ch, r);
% clear r col1 col2 y_upper_lim y_lower_lim t t_start t_end

[zq,rq]=meshgrid(linspace(min(zpos),max(zpos),n),linspace(min(rpos),max(rpos),n));
grid2D=struct('zq',zq,'rq',rq);
clear zq rq


data2D = data2Dcalc(EF, grid2D, n, trange, rpos, zpos, bz, ok);
end
