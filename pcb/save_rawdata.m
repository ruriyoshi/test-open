%1000*128のrawdataをmat形式で保存

%%%offsetなしのrawdataを保存する場合
pathname.rawdata=getenv('rawdata_path'); %保存先
dtacqlist=10005:10007;%IDXではなくdtacqの通し番号

for dtacq=dtacqlist(1,:)
    clear x
    [rawdata]=getvalue_noTF(dtacq);
    save(strcat(pathname.rawdata,'rawdata_noTF_dtacq',num2str(dtacq),'.mat'),'rawdata');
end

%TFoffsetの対応がある場合に、オフセットを引いて保存する場合
% %%%%(1)spread sheetから ログのテーブルを取得してTに格納
% %Github/test-open/getTS6log.mを使用
% DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
% T=getTS6log(DOCID);
% 
% pathname.rawdata=getenv('woTFdata_path'); %保存先
% 
% %%%(2)d_tacqとTFdtacqの両方のデータがある場合に絞り込み
% IDXlist=1;%1:2950;%2452:2950;
% Nan1=isnan(T.d_tacq(IDXlist));
% Nan2=isnan(T.TFdtacq(IDXlist));
% Nan=Nan1|Nan2;
% IDXlist(Nan)=[];
% 
% %%%%(3)指定したshotの解析
% %IDXlist=[2911:2913 2925 2926 2927 2931 2933 2947:2950 2942 2943 2946];
% for IDX=IDXlist(1,:)
%     d_tacq=T.d_tacq(IDX);
%     d_tacqTF=T.TFdtacq(IDX);
%     [rawdata]=getvalue(d_tacq,d_tacqTF);
%     save(strcat(pathname.rawdata,'rawdata_dtacq',num2str(d_tacq),'.mat'),'rawdata');
% end



