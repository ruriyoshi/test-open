function [ts6log]=getTS6log(DOCID)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Google spread sheetをcsvにしてMatlabに変数として読み込�?
% [DOCID]はGoogle spread sheetのリンクのhttps://docs.google.com/spreadsheets/d/のあとの部�?
% [DOCID] A value like '0AmQ013fj5234gSXFAWLK1REgwRW02hsd3c', which is found in your spreadsheet's url: https://docs.google.com/spreadsheets/d/<here>/edit#gid=0.
% 参�??��https://jp.mathworks.com/matlabcentral/fileexchange/39915-getgooglespreadsheet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

loginURL = 'https://www.google.com';
csvURL = ['https://docs.google.com/spreadsheet/ccc?key=' DOCID '&output=csv&pref=2'];
cookieManager = java.net.CookieManager([], java.net.CookiePolicy.ACCEPT_ALL);
java.net.CookieHandler.setDefault(cookieManager);
handler = sun.net.www.protocol.https.Handler;
connection = java.net.URL([],loginURL,handler).openConnection();
connection.getInputStream();

options=weboptions('CharacterEncoding','UTF-8');

filename='ts6log.csv';
outfilename=websave(filename,csvURL,options);
opts=detectImportOptions(filename);
%ts6log = webread(csvURL,options);
ts6log=readtable(outfilename,opts);
end
