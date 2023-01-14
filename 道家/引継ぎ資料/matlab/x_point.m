function [x_point_index] = x_point(B_z,r_probe)

psi_all = zeros(size(B_z)); %psi(6,28,1024)


% get psi
for i = 1:1024
    psi_all(:,:,i) = get_psi(B_z,r_probe,i);
end

max_value_alongz = zeros(1,28,1024); % create array (store max value of psi at each z)
r_direction_index = zeros(1,28,1024); % create array (the index of r-position of max value) 

for i=1:1024
    [M, I] = max(psi_all(:,:,i)); % get the max value 
    max_value_alongz(:,:,i) = M; % max value of psi at each z
    r_direction_index(:,:,i) = I; % the index of r-position of max value
end

x_point_of_z_array = islocalmin(max_value_alongz); % get local min in the form of logical

for i=1:1024
    if nnz(x_point_of_z_array(:,:,i)) == 0  % if there is not local min, set 1 at z_probe(14)
        x_point_of_z_array(1,14,i) = 1;
     end
    if nnz(x_point_of_z_array(:,:,i)) > 1  % if there are some local min, set non minimum index to 0  
        index = find(x_point_of_z_array(:,:,i) == 1);
        disc = zeros(1,length(index));
        for j=1:length(index)
            disc(j) = max_value_alongz(1,index(j),i);
        end
        [Min, I_min] = min(disc);
        all_min_index = index(I_min);
        x_point_of_z_array(:,:,i) = zeros(1,28);
        x_point_of_z_array(:,all_min_index,i) = 1;
    end
end
        

x_point_of_z_index = zeros(1024,1);
for i=1:1024
    x_point_of_z_index(i,1) = find(x_point_of_z_array(:,:,i) == 1);
end
x_point_of_r_index = r_direction_index(x_point_of_z_array);

x_point_index = zeros(size(B_z));
for i=1:1024
    x_point_index(x_point_of_r_index(i),x_point_of_z_index(i),i) = 1;
end

x_point_index = logical(x_point_index);

end