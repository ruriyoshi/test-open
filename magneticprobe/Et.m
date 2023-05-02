function E_t = Et(B_z,r_probe,time)

psi_before = get_psi(B_z,r_probe,time-1);
psi_now = get_psi(B_z,r_probe,time);

E_t = zeros(size(psi_now));
for i = 1:length(r_probe)
    E_t(i,:) = -1/(2*pi*r_probe(i)) * (psi_now(i,:) - psi_before(i,:));
end
E_t = E_t.*1e+6; %単位時間1usによる微分を標準単位系sに直す
end
