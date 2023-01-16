function phase = total_phase(freq)
% get phase from frequency
    phase = zeros(size(freq));
    for k = 1:length(phase)
        phase(k) = sum(freq(1:k))*2*pi;
    end
end