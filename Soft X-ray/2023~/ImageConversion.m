function [] = ImageConversion(foldername)
% path = '/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/SXR_Images/';
% path = strcat(path,foldername);
sourcepath = foldername;% 軟X線画像の保存されているフォルダのパス
if exist(sourcepath,'dir') == 0
    disp('Inadequate path')
    return
end
MyFolderInfo = dir(sourcepath);
outputpath = strcat(sourcepath,'/converted');
if exist(outputpath,'dir') == 0
    mkdir(outputpath);
    NumData = numel(MyFolderInfo);
else
    NumData = numel(MyFolderInfo)-1;
end
MyFolderInfo = MyFolderInfo(4:NumData);
MyFolderDir = MyFolderInfo.folder;
MyFolderName = {MyFolderInfo.name};
FolderInfo = strcat(MyFolderDir,'/',MyFolderName);

for i = 1:numel(FolderInfo)
    if MyFolderInfo(i).isdir
        continue
    end
    figure(1);
    % figure;
    IM = imread(string(FolderInfo(i)));
    imagesc(IM,[50,80]);
    imagesc(IM);
    saveas(gcf,strcat(outputpath,'/',string(MyFolderName(i))));
end
end
