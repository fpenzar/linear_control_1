G_2 = tf([0.13], [0.0005, 0.0115, 0.08, 1]);
bode(G_2)

Ni_vel = 3;
alpha_vel = 0.2;
phase_margin_vel = 60;

[wc, Kp, taui, taud, ok] = findpid(G_2, phase_margin_vel, Ni_vel, alpha_vel);
wc

Cpi = tf([taui, 1], [taui 0]);
Cd = tf([taud, 1], [alpha_vel*taud, 1]);
G_ol = Kp*Cpi*Cd*G_2;
G_cl_fwd = G_ol / (1 + G_ol);

bode(G_ol);
step(G_cl_fwd);

% [Kp_vel, Cpi_vel, Cd_vel, ti_vel, td_vel, w_c_vel, G_ol_vel, G_cl_fwd_vel, G_cl_fdb_vel] = PI_lead(phase_margin_vel, alpha_vel, Ni_vel, G_tilt_vel);
% step(G_cl_fwd_vel)