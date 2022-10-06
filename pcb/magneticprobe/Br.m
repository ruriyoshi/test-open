function B_r = Br(psi,z_space,r_space)
% calculate B_r using B_r = -1/(2*pi*r)*dPsi/dz

B_r = zeros(size(psi));
for i = 1:length(r_space)
    pp = spline(z_space,psi(i,:));
    p_der = fnder(pp,1);
    B_r(i,:) = -1/(2*pi*r_space(i))*ppval(p_der,z_space);
    B_r(i,:) = smooth(z_space,B_r(i,:),5,'moving');
end

end