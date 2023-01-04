direction_list = ['T','R','Z'];

xls_file = 'rms_data_without_resistor50.xls';

col_names = {'RowNames','rms_func_gene', 'rms_amp', 'rms_coil_raw', 'rms_coil_smoothed'};
writecell(col_names,xls_file);

for direction = direction_list

rms_mat = low_n_calib(direction);


%%%%%%%% export to table %%%%%%%%

rowNames = strings(8,1);
for i = 1:8
    formatSpec = "%s-%d";
    rowNames(i,1)= compose(formatSpec, strcat('B',direction), i);
end




export_mat = [rowNames, rms_mat];
writematrix(export_mat,xls_file,'WriteMode','append');


rms_mat = zeros(4,8);
end