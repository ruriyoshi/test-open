function h = plot_psi_multi_pcb_for_axis(data2D,grid2D,fill,fixed_Clayer,shot,n,trange,xpos,max_psi_left_pos, max_psi_right_pos)
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

h = figure('Position', [0 0 800 1500],'visible','off');
%h = figure('Position', [0 0 800 1500],'visible','on');
row = 4;
column = 4;
time_offset = 399;
%time = (482:3:527)-time_offset;
time = (450:3:495)-time_offset;
%time = (420:3:465) - time_offset;
%time = (464:1:473)-time_offset;
%time = 420:4:480;
%time = 450:4:562;
%time = 450:4:470;
%time = (463:1:474)-time_offset;
%time = 450:32:482;
%time = 480;
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
show_xpoint = true;
show_axis_point = true;

psi_mesh_z = grid2D.zq;
psi_mesh_r = grid2D.rq;

% jt

min_color = -1.5*1e+6;
max_color = 1.5*1e+6;
contour_layer_color =  min_color:(max_color-min_color)/50:max_color;


%min_color = -0.03;
%max_color = 0.03;
%contour_layer_color =  min_color:(max_color-min_color)/50:max_color;

if fixed_Clayer
    max_psi = 60e-3;
    min_psi = -20e-3;
    layer_resolution = 0.2e-3;
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
            contourf(psi_mesh_z,psi_mesh_r,color_store(:,:,i),contour_layer,'Fill','off','LineWidth',0.6,'LineColor','white');
        end
        contour(psi_mesh_z,psi_mesh_r,psi_store(:,:,i),contour_layer,'Fill','off','LineWidth',0.6);
    end
    
    if show_probe
        plot(probe_mesh_z,probe_mesh_r,'*','color','k','markersize',1)
    end

    if show_xpoint
        plot(xpos(1,i),xpos(2,i),'bx')%X点
    end
    
    if show_axis_point
        plot(max_psi_left_pos(1,i),max_psi_left_pos(2,i),'r*');
        plot(max_psi_right_pos(1,i),max_psi_right_pos(2,i),'r*');
    end
    

    
    %circle(0.35,0.22,0.0375);
    %circle(-0.35,0.22,0.0375);
    %rectangle('Position',[z_probe(14),r_probe(end),z_probe(15)-z_probe(14),r_probe(1)-r_probe(end)],'EdgeColor','w','FaceColor','w')

    title(strcat(num2str(i+time_offset),' us'),'FontSize',12);
   
    xlabel('z (m)');
    ylabel('r (m)');
    n = n+1;
    %xlim([-0.1 0.1]);
    %ylim([0.15 0.33]);
    daspect([1 1 1])
    hold off
    
end
tiles.TileSpacing = 'tight';
tiles.Padding = 'tight';
cb = colorbar;
cb.Layout.Tile = 'east';
colormap(jet)
axis image
axis tight manual

function h = circle(x,y,r)
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit,"Color",'k');
end
end