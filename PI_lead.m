function [Kp, Cpi, Cd, ti, td, w_c, G_ol, G_cl_fwd, G_cl_fdb] = PI_lead(phase_margin, alpha, Ni, G)
    % Find the new cross-over frequency
    w_c = phaseBalance_equation(deg2rad(phase_margin), alpha, Ni, G);
    % w_c = 15;

    % Find the time constant ti and Cpi
    ti = Ni / w_c;
    Cpi = tf([ti, 1], [ti 0]);

    % Find td and Cd
    td = 1 / (w_c * sqrt(alpha));
    Cd = tf([td, 1], [alpha*td, 1]);

    % Find the Kp
    open_loop = Cpi * Cd * G;
    [mag, ~, ~] = bode(open_loop, w_c); % 'bode' returns magnitude in absolute units (not dB), so no need to convert
    Kp = 1 / squeeze(mag);

    % calculate open loop tf now
    G_ol = Kp*Cpi*Cd*G;

    % calculate closed loop in forward branch
    G_cl_fwd = G_ol / (1 + G_ol);
    % calculate closed loop in backward branch
    G_cl_fdb = (Kp*Cpi*G) / (1 + Kp*Cpi*G*Cd);
end