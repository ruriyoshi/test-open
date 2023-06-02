function h = plot_psi_multi_pcb(data2D,grid2D,fill,fixed_Clayer,shot)
% plot psi in rz plane at multiple times
% input:
%   3d array of double: B_z (r,z,t), offsetted at zero and smoothed
%   1d array of double: r_probe, locations of probes along r
%   1d array of double: z_probe, locations of probes along z
%   3d array of double: B_t (r,z,t), can use B_z to replace it if not used
%   1d array of double: r_probe_Bt, locations of probes along r for Bt
%   1d array of double: z_probe_Bt, locations of probes along z for Bt
%   boolean: fitting, option for enabling fitting and smoothing
%   boolean: fixed_Clayer, option for enabling 
%            true -> fixed value between each contour layer
%            false-> fixed number of total layers = layer_num

% Number of plots in a row and in a column
% It is required that row * column = length(time)

%h = figure('Position', [0 0 1500 1500],'visible','off');
h = figure('Position', [0 0 800 1000],'visible','on');
row = 4;
column = 4;
time_offset = 399;
%time = (482:3:527)-time_offset;
%time = (450:3:495)-time_offset;
time = (457:1:472)-time_offset;
%time = (471:1:480)-time_offset;
%time = (480:3:525)-time_offset;
%time = (420:3:465) - time_offset;
%time = (464:1:473)-time_offset;
%time = 420:4:480;
%time = 450:4:562;
%time = 450:4:470;
%time = 450:32:482;
%time = 480;
%time = 462-time_offset;
%time = 440:3:485;
%time = 440:4:500;
%time = 430:4:490;
%
%time = (480:4:540)-time_offset;
%time = 405:6:495;
%time = 450:10:600;
%time = 490;

tiles = tiledlayout(row,column,'TileSpacing','tight');
layer_num = 30;
show_probe = false;
show_jt = true;
show_Bz = false;
show_Et = false;
show_Bt = false;


psi_mesh_z = grid2D.zq;
psi_mesh_r = grid2D.rq;

% jt

min_color = -3.0*1e+6;
max_color = 3.0*1e+6;
contour_layer_color =  min_color:(max_color-min_color)/50:max_color;


%Bt
%{
min_color = -0.05;
max_color = 0.05;
contour_layer_color =  min_color:(max_color-min_color)/50:max_color;
%}
if fixed_Clayer
    max_psi = 60e-3;
    min_psi = -20e-3;
    layer_resolution = 0.5e-3;
    contour_layer =  min_psi:layer_resolution:max_psi;
else
    contour_layer = layer_num;
end

psi_store = zeros([size(data2D.psi,1:2),length(time)]);
color_store = psi_store;
Br_store = psi_store;
Bz_store = psi_store;

for i = time    
    psi = data2D.psi(:,:,i);
    psi_store(:,:,i) = psi;
    
    if show_jt
        color_store(:,:,i) =  -1.*data2D.Jt(:,:,i);
    end
    if show_Bz
        color_store(:,:,i) = data2D.Bz(:,:,i);
    end
    if show_Et
        color_store(:,:,i) = -1.*data2D.Et(:,:,i);
    end
    if show_Bt
        color_store(:,:,i) = data2D.Bt(:,:,i);
    end
    if not(fill)
        Br_store(:,:,i) = data2D.Br(:,:,i);
        Bz_store(:,:,i) = data2D.Bz(:,:,i);
    end

end


    

j = 1;
for i = time
    nexttile
    hold on
    
    if fill
        contourf(psi_mesh_z,psi_mesh_r,psi_store(:,:,i),contour_layer,'Fill','on');
%         contourf(psi_mesh_z,psi_mesh_r,psi_store(:,:,i),[0 0],'Fill','off','LineWidth',2);
    else
        if (show_jt || show_Et || show_Bz || show_Bt)
            contourf(psi_mesh_z,psi_mesh_r,color_store(:,:,i),contour_layer_color,'LineStyle','none');
            caxis([min_color max_color])
            contour(psi_mesh_z,psi_mesh_r,color_store(:,:,i),[0.7e6 0.64e6 0.53e6],'w-','LineWidth',0.1);
            
        end
        contourf(psi_mesh_z,psi_mesh_r,psi_store(:,:,i),contour_layer,'Fill','off','LineWidth',0.1);
    end
    
    if show_probe
        plot(probe_mesh_z,probe_mesh_r,'*','color','k','markersize',1)
    end



    
    %circle(0.35,0.22,0.0375);
    %circle(-0.35,0.22,0.0375);
    %rectangle('Position',[z_probe(14),r_probe(end),z_probe(15)-z_probe(14),r_probe(1)-r_probe(end)],'EdgeColor','w','FaceColor','w')

    title(strcat(num2str(i+time_offset),' us'),'FontSize',12);
   
    xlabel('z (m)');
    ylabel('r (m)');
    
    %xlim([-0.1 0.1]);
    %ylim([0.1 0.25]);


    daspect([1 1 1])
    j = j+1;
    hold off
    
end
tiles.TileSpacing = 'tight';
tiles.Padding = 'tight';
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Jt (A/m^2)';
%caxis([-8e-3,8e-3])%psi
colormap(jet)
axis image
axis tight manual

function h = circle(x,y,r)
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit,"Color",'k');
end

clear all
end