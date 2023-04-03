function [] = rgw2txt_shot(date_str,shot_str)

% current_folder = strcat('/Users/shinjirotakeda/mountpoint/',date,'/');
folder_directory_rogo = getenv('rogo_path');
current_folder = strcat(folder_directory_rogo,date_str,'/');%strcat('/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/probedata/',date_str,'/');
filename = strcat(current_folder,date_str,shot_str,'.rgw');
rename = strcat(current_folder,date_str,shot_str,'.txt');
copyfile(filename,rename);

end