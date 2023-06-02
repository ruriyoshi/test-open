function phase = cumulative_phase(Ph)
% get cumulative phase from phase
    phase = zeros(size(Ph));
    phase_reference = 0;
    phase(1) = phase_reference + Ph(1);
    for k = 2:length(phase)
        if Ph(k)-Ph(k-1) > pi
            phase_reference = phase_reference - 2*pi;
        elseif Ph(k)-Ph(k-1) < -pi
            phase_reference = phase_reference + 2*pi;
        end
        phase(k) = phase_reference + Ph(k);
    end
end