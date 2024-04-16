function [w_c] = phaseBalance_equation(phase_margin, alpha, Ni, G)
    phase_w_c = rad2deg(phase_margin - asin((1-alpha)/(1+alpha)) - pi + atan2(1, Ni));

    [~, phase, freq] = bode(G);
    phase = squeeze(phase);
    freq = squeeze(freq);
    
    % Calculate the absolute difference between phase_w_c and all phase values
    phase_diff = abs(phase - phase_w_c);
    % Find the index of the minimum difference
    [~, idx] = min(phase_diff);

    % Find the corresponding frequency
    w_c = freq(idx);
end