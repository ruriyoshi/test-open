function [] = plot_psi_xpoint(B_z,r_probe,z_probe,time,fitting,fill,fixed_Clayer,show_probe,show_xo,IDX)
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

figure('Position', [0 0 1500 1500],'visible','on');

column = ceil(sqrt(numel(time))); 
row = ceil(numel(time)/column);

%time = 460:2:490;
%time = 420:10:570;

layer_num = 30;
mesh_z=100;
mesh_r=100;
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
    layer_resolution = 0.0002; %plotする等高線の間隔
    contour_layer =  min_psi_color:layer_resolution:max_psi_color;
else
    contour_layer = layer_num;
end

if fitting
    z_space = linspace(z_probe(1),z_probe(end),mesh_z);
    r_space = linspace(r_probe(1),r_probe(end),mesh_r);
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

if show_xo
    psi_or=zeros(2,length(time));
    psi_oz=psi_or;
    psi_xr=zeros(1,length(time));
    psi_xz=psi_xr;
    [max_psi,max_psi_r]=max(psi_store,[],1); %各時間、各列(z)ごとのpsiの最大値
    max_psi=squeeze(max_psi);
    
    for i = time
        ind=i-time(1)+1;
        max_psi_ind=find(islocalmax(smooth(max_psi(:,i)),'MaxNumExtrema', 2));
        r_ind=max_psi_r(:,:,i);
        if numel(max_psi(min(max_psi_ind),i))==0
            psi_or(1,ind)=NaN;
            psi_oz(1,ind)=NaN;
        else
            o1r=r_ind(min(max_psi_ind));
            psi_or(1,ind)=r_space(o1r);
            psi_oz(1,ind)=z_space(min(max_psi_ind));
        end
        if numel(max_psi(max(max_psi_ind),i))==0
            psi_or(2,ind)=NaN;
            psi_oz(2,ind)=NaN;
        else
            o2r=r_ind(max(max_psi_ind));
            psi_or(2,ind)=r_space(o2r);
            psi_oz(2,ind)=z_space(max(max_psi_ind));
        end
        if numel(find(islocalmin(smooth(max_psi(:,i)),'MaxNumExtrema', 1)))==0
            psi_xr(1,ind)=NaN;
            psi_xz(1,ind)=NaN;
        else
            min_psi_ind=islocalmin(smooth(max_psi(:,i)),'MaxNumExtrema', 1);
            xr=r_ind(min_psi_ind);
            if xr==1 || xr==mesh_r %r両端の場合は検知しない
                psi_xr(1,ind)=NaN;
                psi_xz(1,ind)=NaN;   
            else
                psi_xr(1,ind)=r_space(xr);
                psi_xz(1,ind)=z_space(min_psi_ind);
            end
        end
        if max_psi(1,i)==max(max_psi(1:end/2,i)) %z両端の場合はpsi中心が画面外として検知しない
            psi_or(1,ind)=NaN; %r_space(r_ind(1));
            psi_oz(1,ind)=NaN; %z_space(1);
        end
        if max_psi(end,i)==max(max_psi(end/2:end,i))
            psi_or(2,ind)=NaN; %r_space(r_ind(1));
            psi_oz(2,ind)=NaN; %z_space(end);
        end
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
    
    if show_xo
        ind=i-time(1)+1;
        plot(psi_oz(1,ind),psi_or(1,ind),'ro','markersize',4)
        plot(psi_oz(2,ind),psi_or(2,ind),'ro','markersize',4)
        plot(psi_xz(1,ind),psi_xr(1,ind),'rx','markersize',4)
    end
    
    hold off
%     xlim([-0.1 0.1])
    title(strcat(num2str(i),' us'),'FontSize',12);
    ax = gca;
    
    xlabel('z (m)','FontSize',12);
    ylabel('r (m)','FontSize',12);
    ax.FontSize = 12;
    
    %k = 0.04;
    %ax.TickLength = [k, k]; % Make tick marks longer.
    %ax.LineWidth = 100*k;
    
    n = n+1;
    daspect([1 1 1])
end

sgtitle(strcat('IDX=',num2str(IDX)))

end