function [] = rgw2txt(date)

% Get all rgw files in the current folder
date = num2str(date);
folder_directory_rogo = getenv('rogo_path');
current_folder = strcat(folder_directory_rogo,date,'/');%strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/probedata/',date,'/');
files = dir(strcat(current_folder,'*.rgw'));
% Loop through each
% disp(length(files));
for id = 1:length(files)
    % Get the file name (minus the extension)
    [~,name,~] = fileparts(files(id).name);
    rename = strcat(current_folder,name,'.txt');
    old_name = strcat(current_folder,files(id).name);
    copyfile(old_name,rename);
%     disp(id);
end

end