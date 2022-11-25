function [Bt, Bz] = angle_calibration(Bt_raw, Bz_raw)

alpha = Bt_raw/Bz_raw;
theta = atan(1/alpha);
mat = [cos(theta) sin(theta); -sin(theta) cos(theta)];
B = mat*[Bt_raw; Bz_raw];
Bt = B(1,:);
Bz = B(2,:);

end


