function [] = plot_SXR_multi(B_z,r_probe,z_probe,date,shot,layer,area,start,interval,save,SXRfilename)
% plot SXR emission on psi in rz plane
% input:
%   3d array of double: B_z (r,z,t), offsetted at zero and smoothed
%   1d array of double: r_probe, locations of probes along r
%   1d array of double: z_probe, locations of probes along z
%   integer: date, date of experiment
%   integer: shot, number of shot
%   boolean: area, option for narrowing the reconstruction area
%   boolean: layer, option for changing the contour property
%   integer: start, start time (us)
%   integer: interval, interval time of the framing camera (us)
%   boolean: save, option for saving the reconstruction result
%   string: SXRfilename, name of the SXR image file

filepath = '/Users/shinjirotakeda/Documents/GitHub/SXR_diagnostics/parameters.mat';
N_projection_new = 80;
N_grid_new = 100;
if isfile(filepath)
    load(filepath, 'U1','U2','s1','s2','v1','v2','M','K','range','N_projection','N_grid');
    if N_projection_new ~= N_projection || N_grid_new ~= N_grid
        disp('Different parameters - Start calculation!');
        [U1,U2,s1,s2,v1,v2,M,K,range] = clc_parameters(N_projection_new,N_grid_new);
    end
else
    disp('No parameter - Start calculation!');
    [U1,U2,s1,s2,v1,v2,M,K,range] = clc_parameters(N_projection_new,N_grid_new);
end

times = start:interval:(start+interval*7);
plot_flag = false;

for t = times
    number = (t-start)/interval+1;

    if date <= 210924
        [VectorImage1,VectorImage2] = get_SXRImage(date,number,SXRfilename);
    else
        [VectorImage2,VectorImage1] = get_SXRImage(date,number,SXRfilename);
    end

    EE1 = clc_distribution(M,K,U1,s1,v1,VectorImage1,plot_flag);
    EE2 = clc_distribution(M,K,U2,s2,v2,VectorImage2,plot_flag);

    plot_save_SXR(B_z,r_probe,z_probe,range,date,shot,t,EE1,EE2,area,layer,save);
end

end