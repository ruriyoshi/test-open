function [] = plot_Et_and_psi_at_t(B_z,r_probe,z_probe,t,fitting,fixed_Clayer,show_probe,shot)
% plot psi in rz plane
% input:
%   3d array of double: B_z (r,z,t), offsetted at zero and smoothed
%   1d array of double: r_probe, locations of probes along r
%   1d array of double: z_probe, locations of probes along z
%   integer: t, time of interest (us)
%   boolean: fitting, option for enabling cubic spline fitting 
%   to fit the the 8x28 mesh to 50x50 mesh
%   boolean: fixed_Clayer, option for enabling 
%            true -> fixed value between each contour layer
%            false-> fixed number of total layers = layer_num
%   boolean: show_probe, option for showing location of pickup coils

layer_num = 100;
% a uniform mesh for interpolation
[probe_mesh_z,probe_mesh_r] = meshgrid(z_probe,r_probe);

if fitting
    z_space = linspace(z_probe(1),z_probe(end),50);
    r_space = linspace(r_probe(1),r_probe(end),50);
    [psi_mesh_z,psi_mesh_r] = meshgrid(z_space,r_space);
else
    [psi_mesh_z,psi_mesh_r] = meshgrid(z_probe,r_probe);
end

if fixed_Clayer
    max_psi_color = 0.0400;
    min_psi_color = -0.0200;
    layer_resolution = 0.0002;
    contour_layer =  min_psi_color:layer_resolution:max_psi_color;
else
    contour_layer = layer_num;
end

[Et_mesh_z,Et_mesh_r] = meshgrid(z_probe(1,9:19),r_probe);

figure('Position', [0 0 1200 600],'name',['shot', num2str(shot)])
hold on

psi = get_psi(B_z,r_probe,t);
Et = get_Et(B_z,r_probe);

if fitting
    psi = griddata(z_probe,r_probe,psi,psi_mesh_z,psi_mesh_r,'cubic');
end

% plot psi line
contour(psi_mesh_z,psi_mesh_r,psi,contour_layer,'-k');


% plot Et map
contourf(Et_mesh_z,Et_mesh_r,Et(:,9:19,t),'LineStyle','none');
c = colorbar;
colormap turbo;
%caxis([-200 30]);
c.Label.String = 'Et [V/m]';

g = get(gca);
set(gca,'Children',[g.Children(2),g.Children(1)]);


if show_probe
    plot(probe_mesh_z,probe_mesh_r,'*','color','k','markersize',3);
end
    
title(strcat('Psi and Et :',num2str(t),' us'));
xlabel('z (m)');
ylabel('r (m)');
ax = gca;
ax.FontSize = 18; 
hold off
daspect([1 1 1])
end