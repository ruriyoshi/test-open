function [error_high_value,error_low_value] = TestTomography

NL = false;

% 再構成計算に必要なパラメータを計算するなら読み込む、しない場合も範囲に関しては読み込む
filepath = '/Users/shinjirotakeda/Documents/GitHub/test-open/Soft X-ray/~2022/parameters.mat';
N_projection_new = 80;
N_grid_new = 100;
if isfile(filepath)
    load(filepath, 'gm2d1', 'gm2d2', 'U1', 'U2', 's1', 's2', 'v1', 'v2', 'M', 'K', 'range','N_projection', 'N_grid');
%         load(filepath);
    if N_projection_new ~= N_projection || N_grid_new ~= N_grid
        disp('Different parameters - Start calculation!');
        clc_parameters_old(N_projection_new,N_grid_new);
        load(filepath, 'gm2d1', 'gm2d2', 'U1', 'U2', 's1', 's2', 'v1', 'v2', 'M', 'K', 'range');
    end
else
    disp('No parameter - Start calculation!');
    clc_parameters_old(N_projection_new,N_grid_new);
    load(filepath, 'gm2d1', 'gm2d2', 'U1', 'U2', 's1', 's2', 'v1', 'v2', 'M', 'K', 'range');
end

plot_flag = false;


image_error = 0.0374;
VectorImage = AssumedProfile(gm2d1,N_projection,range,true);
VectorImage_high = VectorImage.*(1+image_error);
VectorImage_low = VectorImage.*(1-image_error);

% whos U1
% whos VectorImage

%         再構成計算
EE = clc_distribution_old(M,K,gm2d1,U1,s1,v1,VectorImage,plot_flag,NL);
EE_high = clc_distribution_old(M,K,gm2d1,U1,s1,v1,VectorImage_high,plot_flag,NL);
EE_low = clc_distribution_old(M,K,gm2d1,U1,s1,v1,VectorImage_low,plot_flag,NL);

range = range./1000;
zmin1 = range(1);
zmax1 = range(2);
rmin = range(5);
rmax = range(6);
r_space_SXR1 = linspace(rmin,rmax,size(EE,1));
z_space_SXR1 = linspace(zmin1,zmax1,size(EE,2));
[SXR_mesh_z1,SXR_mesh_r1] = meshgrid(z_space_SXR1,r_space_SXR1);

figure;contourf(SXR_mesh_z1,SXR_mesh_r1,EE,20);
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
title('original');

figure;contourf(SXR_mesh_z1,SXR_mesh_r1,EE_high,20);
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
title('with positive error');

figure;contourf(SXR_mesh_z1,SXR_mesh_r1,EE_low,20);
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
title('with negative error');

error_high_image = (EE_high - EE)./EE;
error_low_image = (EE_low - EE)./EE;

figure;contourf(SXR_mesh_z1,SXR_mesh_r1,error_high_image,20);
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
title('positive error');

figure;contourf(SXR_mesh_z1,SXR_mesh_r1,error_low_image,20);
c=colorbar;c.Label.String='Intensity [a.u.]';c.FontSize=18;
title('negative error');

N = numel(EE);
EE_max = max(max(EE));
error_high_value = sum(sum((EE_high./EE_max).^2))/N;
error_low_value = sum(sum((EE_low./EE_max).^2))/N;

mean_error_positive = mean(error_high_image,'all')
mean_error_negative = mean(error_low_image,'all')

% plot_save_SXR(B_z,r_probe,z_probe,range,date,shot,t,EE_high,EE_low,show_localmax,show_xpoint,save,filter,NL);

end

function Iwgn = AssumedProfile(gm2d,N_projection,range,plot_flag)

[~, N_g] = size(gm2d);
% N_projection = sqrt(N_p);
N_grid = sqrt(N_g);
m=N_grid;
n=N_grid;
z_0=0;
r_0=-0.3;
z=linspace(-1,1,m);
r=linspace(-1,1,n);

z_grid = linspace(range(2),range(1),m);
r_grid = linspace(range(6),range(5),n);

EE = zeros(m,n);
for i=1:m
    for j=1:n
        r0 = sqrt((z(i)-z_0)^2+(r(j)-r_0)^2);
        r1 = abs(0.5*(z(i)-z_0)+(r(j)-r_0))/sqrt(1.25);
        EE(i,j) = 1*exp(-0.5*r0^2)*exp(-5*r1)+1.5*exp(-r0^2);
    end
end
% size(EE)
EE = EE./max(max(EE));
EE = fliplr(rot90(EE));

%2D matrix is transformed to 1D transversal vector
E = reshape(EE,1,[]);
% whos gm2d
% whos E
I=gm2d*(E)';
Iwgn=awgn(I,10*log10(10),'measured'); % 5 related to 20%; 10 related to 10%;
Iwgn(Iwgn<0)=0;


% 1D column vector is transformed to 2D matrix
n_p = N_projection;
II = zeros(n_p);
IIwgn = zeros(n_p);
k=FindCircle(n_p/2);

II(k) = I;
IIwgn(k) = Iwgn;

% I_check = II(75,:);
% j = 1:numel(I_check);
% figure;plot(j,I_check);

if plot_flag
    figure;imagesc(z_grid,r_grid,EE);c=colorbar('Ticks',[0.1,0.5,1]);
    c.Label.String='Assumed Intensity [a.u.]';xlabel('Z [mm]');ylabel('R [mm]');
    axis xy
    ax = gca;
    ax.XDir = 'reverse';
    figure;imagesc(II);c=colorbar('Ticks',[0,20,40]);
    c.Label.String='Assumed Intensity [a.u.]';xlabel('Z Pixels');ylabel('R Pixels');
    figure;imagesc(IIwgn);c=colorbar('Ticks',[0,20,40]);
    c.Label.String='Assumed Intensity [a.u.]';xlabel('Z Pixels');ylabel('R Pixels');
end

Iwgn = Iwgn.';

end

function k = FindCircle(L)
R = zeros(2*L);
for i = 1:2*L
    for j = 1:2*L
        R(i,j) = sqrt((L-i+0.5)^2+(j-L-0.5)^2);
    end
end
% figure;imagesc(R)
k = find(R<L);
end
