function [] = ImageConversion(foldername)
% path = '/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/SXR_Images/';
% path = strcat(path,foldername);
% お試し編集
path = foldername;% 軟X線画像の保存されているフォルダのパス
if exist(path,'dir') == 0
    disp('Inadequate path')
    return
end
MyFolderInfo = dir(path);
path2 = strcat(path,'/converted');
if exist(path2,'dir') == 0
    mkdir(path2);
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
    saveas(gcf,strcat(path2,'/',string(MyFolderName(i))));
end
end
