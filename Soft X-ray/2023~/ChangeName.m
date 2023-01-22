date = 211224;
rootDir = strcat('/Users/shinjirotakeda/mountpoint/',num2str(date));% rgwファイルが存在するフォルダのパス
if(rootDir==0)
    disp('処理を終了します');
    return
end

origDir = pwd;% カレントディレクトリをバッファ
cd(rootDir);% カレントディレクトリを変更
fileList = dir('*.rgw');% カレントディレクトリとサブディレクトリを検索
ShotTimeList = cell2mat({fileList.datenum});
cd(origDir);% カレントディレクトリを戻す

% path = strcat("/Users/shinjirotakeda/SXR",num2str(date));% 軟X線画像が保存されているフォルダのパス
path = strcat('/Volumes/experiment/results/X-ray/',num2str(date));% 軟X線画像が保存されているフォルダのパス
Info = dir(path);
Dates = {Info.datenum};
DatesNum = transpose(cell2mat(Dates));


for i = 1:numel(ShotTimeList)
    Idx = knnsearch(DatesNum,ShotTimeList(i));
    dis = abs(DatesNum(Idx)-ShotTimeList(i));
    if dis > 0.0015
        continue;
    end
    Source = strcat(path,'/',Info(Idx).name);
    Destination = strcat(path,'/shots/',strtok(fileList(i).name,'.'),'.tif');
    copyfile(Source,Destination);
end