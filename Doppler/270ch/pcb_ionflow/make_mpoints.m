%�v���_�z��𐶐�
function [mpoints] = make_mpoints(n_CH,min_r,int_r,n_z,min_z,int_z)
mpoints.n_CH = n_CH;
mpoints.n_r = n_CH/4;
mpoints.n_z = n_z;
mpoints.r = zeros(mpoints.n_r,mpoints.n_z);%�x�N�g���v���b�gr���W1���data1�A2���data2
mpoints.z = zeros(mpoints.n_r,mpoints.n_z);%�x�N�g���v���b�gz���W1���data1�A2���data2
for i = 1:mpoints.n_z
    mpoints.r(:,i) = linspace(min_r,min_r+(mpoints.n_r-1)*int_r,mpoints.n_r);%datai��r���W
    mpoints.z(:,i) = min_z+int_z*(i-1);%datai��z���W
end
