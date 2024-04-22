% Assuming you have a transfer function 'sys'
[mag, phase, w] = bode(G_tilt_vel);
phase = phase - 360;

% Unwrap phase data to avoid jumps of 360 degrees
% phase_unwrapped = unwrap(phase);

% Plot the Bode plot manually
figure;
subplot(2,1,1);
semilogx(w, 20*log10(squeeze(mag)));
title('Magnitude (dB)');

subplot(2,1,2);
semilogx(w, squeeze(phase)); % Convert to degrees
title('Phase (deg)');