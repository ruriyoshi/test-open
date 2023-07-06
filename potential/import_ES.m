%% テキスト ファイルからのデータのインポート
% 次のテキスト ファイルからデータをインポートするスクリプト:
%
%    ファイル名: X:\results\ts-3u\210429\ES_No.1.csv
%
% MATLAB からの自動生成日: 2021/07/19 07:39:54

%% インポート オプションの設定およびデータのインポート
function ESdata= import_ES(filename)
opts = delimitedTextImportOptions("NumVariables", 22);

% 範囲と区切り記号の指定
opts.DataLines = [4600, 5000];
opts.Delimiter = ",";

% 列名と型の指定
opts.VariableNames = ["timeus", "ch1", "ch2", "ch3", "ch4", "ch5", "ch6", "ch7", "ch8", "ch9", "ch10", "ch11", "ch12", "ch13", "ch14", "ch15", "ch16", "ch17", "ch18", "ch19", "ch20", "ch21"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% ファイル レベルのプロパティを指定
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% データのインポート
ESNo = readtable(filename, opts);

%% 出力型への変換
ESdata = table2array(ESNo);
ESdata = downsample(ESdata ,10) ;
%% 一時変数のクリア
clear opts
end