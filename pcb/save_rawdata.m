%1000*128のrawdataをmat形式で保存

%%%%(1)spread sheetから ログのテーブルを取得してTに格納
%Github/test-open/getTS6log.mを使用
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);

pathname.rawdata='/Users/mgar/rawdata_a038/'; %保存先

IDXlist=2452:2950;
% Nan_id=isnan(T.d_tacq(IDXlist));
Nan1=isnan(T.d_tacq(IDXlist));
Nan2=isnan(T.TFdtacq(IDXlist));
Nan=Nan1|Nan2;
IDXlist(Nan)=[];


%%%%(3)指定したshotの解析
%IDXlist=[2911:2913 2925 2926 2927 2931 2933 2947:2950 2942 2943 2946];
for IDX=IDXlist(1,:)
    d_tacq=T.d_tacq(IDX);
    d_tacqTF=T.TFdtacq(IDX);
    [rawdata]=getvalue(d_tacq,d_tacqTF);
    save(strcat(pathname.rawdata,'rawdata_dtacq',num2str(d_tacq),'.mat'),'rawdata');
end