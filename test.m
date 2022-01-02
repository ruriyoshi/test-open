clear
%%%%%%%%%%%%%%%%%%%%%%%%
%%　参考　https://jp.mathworks.com/help/matlab/matlab_prog/access-data-in-a-table.html
%%%%%%%%%%%%%%%%%%%%%%%%
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';
T=getTS6log(DOCID);% ログのテーブルを取得

 node='date';  % 検索する列
 pat=211224;   % 検索パターン（数値なら一致検索、文字なら含む検索）

sortedT=searchlog(T,node,pat)% ログのテーブルから当てはまるものを抽出した新しいテーブルを作成



function output=searchlog(inputlogtable,node,pat)
%%名前,commentで検索(~を含むの検索)
if ischar(pat)
rows = contains( inputlogtable.(node),pat,'IgnoreCase',true); 
%%日付,shotで検索(一致)
else
 rows = inputlogtable.(node) == pat ;
end
output=inputlogtable(rows,:);
end

