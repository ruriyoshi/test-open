function [j_t,z_space,r_space] = jt(B_z,z_probe,r_probe,time)
% calculate toroidal current J_t
% input:
%   3d array of double: B_z (r,z,t)
%   1d array of double: z_probe, locations of probes along z
%   1d array of double: r_probe, locations of probes along r
%   integer: time, time of interest

mu0   = 4*pi*10^(-7);   % vacuum permeability (H/m)

z_space = linspace(z_probe(1),z_probe(end),50);
r_space = linspace(r_probe(1),r_probe(end),50);
[psi_mesh_z,psi_mesh_r] = meshgrid(z_space,r_space);

psi = get_psi(B_z,r_probe,time);

B_z_smoothed = B_z(:,:,time);
for j = 1:size(B_z,1)
    B_z_smoothed(j,:) = smooth(z_probe,B_z(j,:,time),3,'moving');
end

psi = griddata(z_probe,r_probe,psi,psi_mesh_z,psi_mesh_r,'cubic');
B_z = griddata(z_probe,r_probe,B_z_smoothed,psi_mesh_z,psi_mesh_r,'cubic');
B_r = Br(psi,z_space,r_space);

B_z_der_in_r = zeros(size(psi));
B_r_der_in_z = zeros(size(psi));

for i = 1:length(r_space)
    B_r_temp = smooth(z_space,B_r(i,:));
    pp = spline(z_space,B_r_temp);
    p_der = fnder(pp,1);
    B_r_der_in_z(i,:) = ppval(p_der,z_space);
    B_r_der_in_z(i,:) = smooth(z_space,B_r_der_in_z(i,:));
end

for i = 1:length(z_space)
    B_z_temp = smooth(r_space,B_z(:,i));
    pp = spline(r_space,B_z_temp);
    p_der = fnder(pp,1);
    B_z_der_in_r(:,i) = ppval(p_der,r_space);
    B_r_der_in_z(:,i) = smooth(r_space,B_r_der_in_z(:,i));
end

j_t = 1/mu0 * (B_r_der_in_z - B_z_der_in_r);
%{
figure('Position', [0 0 1500 1500],'visible','on');
subplot(3,2,1)
contourf(psi_mesh_z,psi_mesh_r,psi,50);
title('Psi')
daspect([1 1 1])
subplot(3,2,2)
contourf(psi_mesh_z,psi_mesh_r,B_z,50);
title('B_z')
daspect([1 1 1])
subplot(3,2,3)
contourf(psi_mesh_z,psi_mesh_r,B_z_der_in_r,50);
title('dBz/dr')
daspect([1 1 1])
subplot(3,2,4)
contourf(psi_mesh_z,psi_mesh_r,B_r,50);
title('B_r')
daspect([1 1 1])
subplot(3,2,5)
contourf(psi_mesh_z,psi_mesh_r,B_r_der_in_z,50);
title('dBr/dz')
daspect([1 1 1])
subplot(3,2,6)
contourf(psi_mesh_z,psi_mesh_r,j_t,50);
title('j_t')
daspect([1 1 1])
%}
end