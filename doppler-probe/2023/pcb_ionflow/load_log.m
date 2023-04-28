%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%実験日(date)を含む実験ログを保存して
%実験日のログ位置(行番号)を取得
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [exp_log,begin_row,end_row] = load_log(date)
%時短用(固定値で良い)
load_s = 5300;%実験ログ読み始め行番号(230309~)
load_f = 10000;%実験ログ読み終わり行番号(固定)
%実験ログがcdになければ最新の実験ログをwebから取得
if exist('exp_log.xlsx','file')
    FileInfo = dir('exp_log.xlsx');
    DateNum = FileInfo.datenum;
    FileDate = datetime(DateNum, 'ConvertFrom', 'datenum', 'Format', 'yyMMdd');
    if str2double(string(FileDate)) > date
        disp(append(string(FileDate), '更新の実験ログを使用します。'));
    else
        save_log();
        disp('最新の実験ログ(exp_log.xlsx)を保存しました。');
    end
else
    save_log();
    disp('最新の実験ログ(exp_log.xlsx)を保存しました。');
end
%実験ログ中の実験日に対応する範囲を特定
exp_log = readmatrix('exp_log.xlsx','Sheet','log','Range', ['A' num2str(load_s) ':AR' num2str(load_f)]);
[n_row,~] = size(exp_log);
begin_row = find(exp_log(:,3) == date,1,"first");%実験日の最初のshotの行番号を取得
if isempty(begin_row)
    error('実験日が実験ログ中に存在しません。')
end
end_row = begin_row;
while end_row<n_row && exp_log(end_row+1,4) && isnan(exp_log(end_row+1,3))... 
        || (exp_log(end_row+1,3) == date)%日付がNaN or date&&shot番号が記入済=実験日のshot
    end_row = end_row+1;
end%実験日の最後のshotの行番号を取得
end
