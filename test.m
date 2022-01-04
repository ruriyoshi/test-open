clear
%%%%%%%%%%%%%%%%%%%%%%%%
%%　参考　https://jp.mathworks.com/help/matlab/matlab_prog/access-data-in-a-table.html
%%%%%%%%%%%%%%%%%%%%%%%%


DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';
T=getTS6log(DOCID);% ログのテーブルを取得

 node='date';  % 検索する列の名前. T.Properties.VariableNamesで一覧表示できる
 pat=211224;   % 検索パターン（数値なら一致検索、文字なら含む検索）　

searchlog(T,node,pat)% ログのテーブルから当てはまるものを抽出した新しいテーブルを作成

%logからshot,TFshot,などのオペレーション設定を読む
prompt = 'number: ';%表示したいショットのnumberの値（通し番号）を入力
IDX=input(prompt) ;
date=T.date(IDX);
shot=T.shot(IDX);
TF_shot = T.TFoffset(IDX) ;
offset_TF = true;

%%ファイルへのパスを作る
%それぞれのPCから共有フォルダまでのパスはそれぞれ異なるので各自で設定
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス
pathname.fourier='I:';%md0までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath
%共有フォルダ以下から目的ショットのファイルを探す
% filepath.mrd=strcat(pathname.ts3u, '\', string(date),'\' ...
%     ,string(date),num2str(shot,'%03i'),'.mrd');
% filepath.tfmrd=strcat(pathname.ts3u, '\', string(date),'\' ...
%     ,string(date),num2str(TF_shot,'%03i'),'.mrd');
filepath.rgw=strcat(pathname.ts3u, '\', string(date),'\' ...
    ,string(date),num2str(shot,'%03i'),'.rgw');
filepath.D288=dir(strcat(pathname.fourier,'\Doppler\288CH\20',string(date),'\*shot',num2str(shot),'*.asc'));
filepath.Dhighspeed=dir(strcat(pathname.NIFS,'\Doppler\Photron\',string(date),'\**\*shot',num2str(shot),'*.tif'));
filepath.SXR=strcat(pathname.NIFS,'\X-ray\',string(date),'\shots\',string(date),num2str(shot,'%03i'),'.tif');

cd magneticprobe
[B_z,r_probe,z_probe,ch_dist,B_z_return,data_return,shot_num] = get_B_z(date,TF_shot,shot,offset_TF,T.EF_A_(IDX),pathname.ts3u);
B_z = B_z([2,3,4,6,7,8],2:end,:);
data = data([2,3,4,6,7,8],2:end,:);
z_probe = z_probe(2:end);
ch_dist = ch_dist([2,3,4,6,7,8],2:end);
r_probe = r_probe([2,3,4,6,7,8]);





function output=searchlog(inputtable,node,pat)
%%名前,commentで検索(~を含むの検索)
if ischar(pat)
rows = contains( inputtable.(node),pat,'IgnoreCase',true); 
%%日付,shotなどで検索(一致検索)
else
 rows = inputtable.(node) == pat ;
end
output=inputtable(rows,:);
end

