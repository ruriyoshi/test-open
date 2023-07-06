function [] = plot_pcb_SXR_multi(psi,rq,zq,date,shot,layer,area,start,interval,save,SXRfilename)
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

if date == 210309
    if shot <= 15
        start = 440;
%         exposure = 5;
        interval = 10;
    else
        start = 450;
%         exposure = 2;
        interval = 5;
    end
elseif date == 210310
    start = 450;
%     exposure = 2;
    interval = 5;
elseif date == 210319
    if shot == 4
        start = 450;
%         exposure = 2;
        interval = 5;
    elseif shot <= 6
        start = 450;
%         exposure = 5;
        interval = 10;
    else
        start = 425;
%         exposure = 5;
        interval = 10;
    end
end

% if exposure == 2
%     interval = 5;
% elseif exposure == 5
%     interval = 10;
% end

times = start:interval:(start+interval*7);

for t = times
    plot_pcb_SXR_at_t(psi,rq,zq,date,shot,t,layer,area,start,interval,save,SXRfilename)%%ここは更新必要
end
end


