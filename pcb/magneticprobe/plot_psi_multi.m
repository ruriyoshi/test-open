function [] = plot_psi_multi(B_z,r_probe,z_probe,time,fitting,fill,fixed_Clayer,show_probe)

%磁気面を指定した時間の分だけ 描く関数
% plot psi in rz plane at multiple times
% input:
%   3d array of double: B_z (r,z,t), offsetted at zero and smoothed
%   1d array of double: r_probe, locations of probes along r
%   1d array of double: z_probe, locations of probes along z
%   boolean: fitting, option for enabling cubic spline fitting
%   boolean: fixed_Clayer, option for enabling 
%            true -> fixed value between each contour layer
%            false-> fixed number of total layers = layer_num
%   boolean: show_probe, option for showing location of pickup coils

% Number of plots in a row and in a column
% It is required that row * column = length(time)

%data画像保存先
pathname.save='C:\Users\uswk0\OneDrive\デスクトップ\data\movie_out\'; %保存先


figure('Position', [0 0 1500 1500],'visible','on');

column = ceil(sqrt(numel(time))); 
row = ceil(numel(time)/column);

%time = 460:2:490;
%time = 420:10:570;

layer_num = 30;
%show_probe = true;
show_jt = true;
show_Bz = false;
show_Et = false;
if (show_jt || show_Et || show_Bz)
    fill=false;
end

% mesh for showing probe
[probe_mesh_z,probe_mesh_r] = meshgrid(z_probe,r_probe);

if fixed_Clayer
    max_psi_color = 0.0400;
    min_psi_color = -0.0200;
    layer_resolution = 0.0002;
    contour_layer =  min_psi_color:layer_resolution:max_psi_color;
else
    contour_layer = layer_num;
end

if fitting
    z_space = linspace(z_probe(1),z_probe(end),50);
    r_space = linspace(r_probe(1),r_probe(end),50);
    [psi_mesh_z,psi_mesh_r] = meshgrid(z_space,r_space);
else
    z_space = z_probe;
    r_space = r_probe;
    [psi_mesh_z,psi_mesh_r] = meshgrid(z_space,r_space);
end

psi_store = zeros(length(r_space),length(z_space),length(time));
color_store = psi_store;
Br_store = psi_store;
Bz_store = psi_store;

for i = time    
    if not(fitting)
        psi = get_psi(B_z,r_probe,i);
        psi_store(:,:,i) = psi;
    else
        for j = 1:size(B_z,1)
            psi(j,:) = smooth(z_probe,B_z(j,:,i),'lowess');
        end
        for j = 1:size(B_z,2)
            psi(:,j) = smooth(r_probe,B_z(:,j,i),'lowess');
        end
        
        psi = get_psi(B_z,r_probe,i);
        
        psi_store(:,:,i) = griddata(z_probe,r_probe,psi,psi_mesh_z,psi_mesh_r,'cubic');
        
        
        for j = 1:size(psi_store(:,:,i),1)
            psi_store(j,:,i) = smooth(z_space,psi_store(j,:,i),'lowess');
        end
        for j = 1:size(psi_store(:,:,i),2)
            psi_store(:,j,i) = smooth(z_space,psi_store(:,j,i),'lowess');
        end
        
    end
    
    if show_jt
        [j_t,z_space_jt,r_space_jt] = jt(B_z,z_probe,r_probe,i);
        color_store(:,:,i) = griddata(z_space_jt,r_space_jt,j_t,psi_mesh_z,psi_mesh_r,'v4');
    end
    
    if show_Bz
        Bz = B_z(:,:,i);
        color_store(:,:,i) = griddata(z_probe,r_probe,B_z(:,:,i),psi_mesh_z,psi_mesh_r,'v4');
    end
    
    if show_Et
        E_t = Et(B_z,r_probe,i);
        color_store(:,:,i) = griddata(z_probe,r_probe,E_t,psi_mesh_z,psi_mesh_r,'v4');
    end
    
    if not(fill)
        Br_store(:,:,i) = Br(psi_store(:,:,i),z_space,r_space);
        Bz_store(:,:,i) = griddata(z_probe,r_probe,B_z(:,:,i),psi_mesh_z,psi_mesh_r,'v4');
    end
        
end
%{
if fixed_Clayer
    minx = min(psi_store,[],'all');
    maxx = max(psi_store,[],'all');
    contour_layer =  minx:(maxx-minx)/40:maxx;
    %contour_layer =  minx:0.001:maxx;
else
    contour_layer = layer_num;
end
%}

n = 1;
min_color = min(color_store(:));
max_color = max(color_store(:));
contour_layer_color =  min_color:(max_color-min_color)/50:max_color;
for i = time
    
    subplot(row,column,n);
    hold on

    
   
    
    if fill
        contourf(psi_mesh_z,psi_mesh_r,psi_store(:,:,i),contour_layer,'Fill','on');
    else
        if (show_jt || show_Et || show_Bz)
            contourf(psi_mesh_z,psi_mesh_r,color_store(:,:,i),contour_layer_color,'LineStyle','none');
        end
        contourf(psi_mesh_z,psi_mesh_r,psi_store(:,:,i),contour_layer,'Fill','off','LineWidth',1);
    end
    
    if show_probe
        plot(probe_mesh_z,probe_mesh_r,'*','color','k','markersize',1)
    end
    
    hold off
    
    title(strcat(num2str(i),' us'),'FontSize',12);
    ax = gca;
    if i == 480
    xlabel('z (m)','FontSize',12);
    ylabel('r (m)','FontSize',12);
    ax.FontSize = 12;
    else
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    end


    %movieの作成
    %saveas(gcf,strcat(pathname.save,'\time',num2str(i),'_shot',num2str(date),'.jpg'))
    %close
    
    %k = 0.04;
    %ax.TickLength = [k, k]; % Make tick marks longer.
    %ax.LineWidth = 100*k;
    
    n = n+1;
    daspect([1 1 1])
end
end