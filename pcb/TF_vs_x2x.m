%%%%%%%%%%%%%%%%%%%%%%%%
%TF vs XtoXの解析
%%%%%%%%%%%%%%%%%%%%%%%%

%%ファイルへのパスを作る
%それぞれのPCから共有フォルダまでのパスはそれぞれ異なるので各自で設定

%%%%(1)spread sheetから ログのテーブルを取得してTに格納
%Github/test-open/getTS6log.mを使用
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);

%%%%%ここが各PCのパス
%環境変数を設定していない場合はパスを''内に全て記入する（使用しないパスは空白''で良い）
% pathname.ts3u='ts3u_path';%old-koalaのts-3uまでのパス（mrdなど）
% pathname.fourier='fourier_path';%fourierのmd0（データックのショットが入ってる）までのpath
% pathname.NIFS='NIFS_path';%resultsまでのpath（ドップラー、SXR）
% pathname.save='/Users/mgar/pcb_save'; %保存先

%%%%(2)ログから解析したいデータを検索
%Github/test-open/searchlog.mを使用

% node='date';  % 【input】検索する列の名前. T.Properties.VariableNamesで一覧表示できる
%  pat=211223;   % 【input】検索パターン（数値なら一致検索、文字なら含む検索）　
% 
% searchlog(T,node,pat); % ログのテーブルから当てはまるものを抽出した新しいテーブルを作成

%%%%(3)指定したshotの解析
IDXlist=2870:2920; %【input】テーブルから解析したいshot番号を抽出して入力

figure
for i=1:size(IDXlist,2)
    IDX=IDXlist(1,i);
    [x2x,v_TF,~]=cal_x2x(T,IDX);
    if x2x>0
        plot(v_TF,x2x,"ro");
    end
    clear v_TF x2x
    hold on
end
hold off
title(strcat('IDX=',num2str(IDXlist(1,1)),':',num2str(IDXlist(1,i))))
xlabel('TF [kV]')
ylabel('X to X [m]')

%%%%%%%%%%%%%%%%%%%%%%%%
%以下、local関数(getinput, cal_x2x)
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%getinput:実験ログ（T）からshot番号（IDX）におけるオペレーションの値を出力
%出力は構造体に変更しても良いかも
function [date, shot, TF_shot, offset_TF, i_EF, start, v_TF, d_tacq, d_tacqTF, trange, n] = getinput(T,IDX)
date=T.date(IDX);
shot=T.shot(IDX);
TF_shot=T.TFoffset(IDX);
offset_TF=isfinite(TF_shot);

if isnan(T.EF_A_(IDX))%%NaNでないことを確認（ログが空白だとNaNになる）
    i_EF=150;
else  %NaNなら150をとりあえず代入、記入されているときはその値を使う
    i_EF=T.EF_A_(IDX);
end

start=T.Period_StartTime_(IDX);
v_TF=T.TF_kV_(IDX); %TF[kV]

d_tacq=T.d_tacq(IDX);
d_tacqTF=T.TFdtacq(IDX);

%trange=460:490;
trange=460:480;

n=50; %rz方向のメッシュ数
end

%%%cal_x2x:X点間の距離計算
function [x2x,v_TF,t_max]=cal_x2x(T,IDX)
[date, ~, ~, ~, i_EF, ~, v_TF, d_tacq, d_tacqTF, trange, n] = getinput(T,IDX);%Tのテーブルから入力のリストを出力

[grid2D, data2D] = pcbdata(date, d_tacq, d_tacqTF, trange, [], n,i_EF);

if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end
   
%%%midplaneとかO点、X点を探す
[psimid,~]=min(data2D.psi,[],2);
%[opoint,p]=islocalmin(psimid,1);
[xpoint,~]=islocalmax(psimid,1);
%[xp_psi,maxxp]=max(squeeze(psimid),[],1);
% onum=squeeze(sum(opoint,1));
% trange(onum~=0)
    %maxrange=2e6;

%%%460-490usの間でX点が2つかつ間隔が最大のものを探す
xnum=zeros(1,size(trange,2));
x2x=0;
t_max=0;
for i=1:size(trange,2)
    xnum(1,i)=sum(xpoint(:,1,i));
    if xnum(1,i)==2
        rq=grid2D.rq(:,1);
        x_rq=rq(xpoint(:,1,i));
        if max(x_rq)-min(x_rq)>x2x
            x2x=max(x_rq)-min(x_rq);
            t_max=459+i;
        end
    end  
end

end

