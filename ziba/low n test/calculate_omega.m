function Omega = calculate_omega(Ph)
    Omega = zeros(size(Ph));
    for j = 1:length(Ph)-1
        Omega(j) = Ph(j+1)-Ph(j);
    end
    Omega(end) = Omega(end-1);
end