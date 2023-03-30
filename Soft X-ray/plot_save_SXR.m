function plot_save_SXR(B_z,r_probe,z_probe,range,date,shot,t,EE1,EE2,area,layer,save)
% plot SXR emission and psi in rz plane and save the result
% input:
%   3d array of double: B_z (r,z,t), offsetted at zero and smoothed
%   1d array of double: r_probe, locations of probes along r
%   1d array of double: z_probe, locations of probes along z
%   1d array of integer: range, [zmin1,zmax1,zmin2,zmax2,rmin,rmax]
%   integer: date, date of the experiment
%   integer: shot, shot number
%   integer: t, time of interest (us)
%   2d array of double: EE1,EE2, local distribution of SXR emission
%   boolean: area, option for narrowing the plot area
%   boolean: layer, option for changing the contour property
%   boolean: save, option for saving the reconstruction result

range = range./1000;
% [zmin1,zmax1,zmin2,zmax2,rmin,rmax] = range;
zmin1 = range(1);
zmax1 = range(2);
zmin2 = range(3);
zmax2 = range(4);
rmin = range(5);
rmax = range(6);
r_space_SXR = linspace(rmin,rmax,size(EE1,1));
z_space_SXR1 = linspace(zmin1,zmax1,size(EE1,2));
z_space_SXR2 = linspace(zmin2,zmax2,size(EE2,2));

if area %この場合は範囲を制限する
    r_range = find(r_probe(end)<=r_space_SXR & r_space_SXR<=r_probe(1));   
    z_range1 = find(-0.12<=z_space_SXR1 & z_space_SXR1<=0.20);
    z_range2 = find(-0.20<=z_space_SXR2 & z_space_SXR2<=0.12);
else %この場合は制限しないor広めに制限
    r_range = find(0.055<=r_space_SXR & r_space_SXR<=0.3);
    z_range1 = find(-0.20<=z_space_SXR1 & z_space_SXR1<=0.20);
    z_range2 = find(-0.20<=z_space_SXR2 & z_space_SXR2<=0.20);
%     r_range = 1:numel(r_space_SXR);
%     z_range1 = 1:numel(z_space_SXR1);
%     z_range2 = 1:numel(z_space_SXR2);
end

r_space_SXR = r_space_SXR(r_range);
z_space_SXR1 = z_space_SXR1(z_range1);
z_space_SXR2 = z_space_SXR2(z_range2);

EE1 = EE1(r_range,z_range1);
EE2 = EE2(r_range,z_range2);

z_space = linspace(z_probe(1),z_probe(end),50);
r_space = linspace(r_probe(1),r_probe(end),50);
[psi_mesh_z,psi_mesh_r] = meshgrid(z_space,r_space);

psi = get_psi(B_z,r_probe,t+2);
psi = griddata(z_probe,r_probe,psi,psi_mesh_z,psi_mesh_r,'cubic');
if layer
    psi_min = min(min(psi));
    psi_max = max(max(psi));
    contour_layer = linspace(psi_min,psi_max,20);
else
    max_psi_color = 0.0400;
    min_psi_color = -0.0200;
    layer_resolution = 0.0002;
    contour_layer =  min_psi_color:layer_resolution:max_psi_color;
end

f = figure;
f.Units = 'normalized';
f.Position = [0.1,0.2,0.8,0.4];
pos1 = [0.07,0.2,0.35,0.6];
pos2 = [0.58,0.2,0.35,0.6];

subplot('Position',pos1);
[SXR_mesh_z1,SXR_mesh_r1] = meshgrid(z_space_SXR1,r_space_SXR);
[~,h1] = contourf(SXR_mesh_z1,SXR_mesh_r1,EE1,20);
h1.LineStyle = 'none';
caxis([0,0.2]);
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
hold on
[~,hp1]=contourf(psi_mesh_z,psi_mesh_r,psi,contour_layer,'-k','Fill','off');
if layer
    hp1.LineWidth = 1.5;
end
title(strcat(num2str(t),' us'));
xlabel('z [m]');
ylabel('r [m]');
ax = gca;
ax.FontSize = 18; 
hold off

subplot('Position',pos2);
[SXR_mesh_z2,SXR_mesh_r2] = meshgrid(z_space_SXR2,r_space_SXR);
[~,h2] = contourf(SXR_mesh_z2,SXR_mesh_r2,EE2,20);
h2.LineStyle = 'none';
% caxis([0,0.15]);
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
hold on
[~,hp2]=contourf(psi_mesh_z,psi_mesh_r,psi,contour_layer,'-k','Fill','off');
if layer
    hp2.LineWidth = 1.5;
end
title(strcat(num2str(t),' us'));
xlabel('z [m]');
ylabel('r [m]');
ax = gca;
ax.FontSize = 18; 
hold off

if save
    pathname = '/Users/shinjirotakeda/OneDrive - The University of Tokyo/Documents/ReconstructionResults/';
    if area
        foldername = strcat(pathname,num2str(date),'/shot',num2str(shot));
    else
        foldername = strcat(pathname,num2str(date),'/shot',num2str(shot),'_wide');
    end
    if exist(foldername,'dir') == 0
        mkdir(foldername);
    end
    filename = strcat('/shot',num2str(shot),'_',num2str(t),'us.png');
    saveas(gcf,strcat(foldername,filename));
end
end