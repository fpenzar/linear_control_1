%% Simscape multibody model og Regbot in balance
% initial setup with motor velocity controller 
% this is intended as simulation base for balance control.
%
close all
clear

%% Simulink model name
model='regbot_1mg';

%% parameters for REGBOT
% motor
RA = 3.3/2;    % ohm (2 motors)
JA = 1.3e-6*2; % motor inertia
LA = 6.6e-3/2; % rotor inductor (2 motors)
BA = 3e-6*2;   % rotor friction
Kemf = 0.0105; % motor constant
Km = Kemf;
% køretøj
NG = 9.69; % gear
WR = 0.03; % wheel radius
Bw = 0.155; % wheel distance
% 
% model parts used in Simulink
mmotor = 0.193;   % total mass of motor and gear [kg]
mframe = 0.32;    % total mass of frame and base print [kg]
mtopextra = 0.97 - mframe - mmotor; % extra mass on top (charger and battery) [kg]
mpdist =  0.10;   % distance to lit [m]
% disturbance position (Z)
pushDist = 0.1; % relative to motor axle [m]

%% Load the variables
% HOMEMADE
load('control_variables.mat');

%% wheel velocity controller (no balance) PI-regulator
% sample (usable) controller values
Kpwv = 15;     % Kp
tiwv = 0.05;   % Tau_i
Kffwv = 0;     % feed forward constant
startAngle = 10;  % tilt in degrees at time zero
twvlp = 0.005;    % velocity noise low pass filter time constant (recommended)

%% Estimate transfer function for base system using LINEARIZE
% Motor volatge to wheel velocity (wv)
load_system(model);
open_system(model);
% define points in model
ios(1) = linio(strcat(model,'/Limit9v'),1,'openinput');
ios(2) = linio(strcat(model, '/wheel_vel_filter'),1,'openoutput');
% attach to model
setlinio(model,ios);
% Use the snapshot time(s) 0 seconds
op = [0];
% Linearize the model
sys = linearize(model,ios,op);
% get transfer function
[num,den] = ss2tf(sys.A, sys.B, sys.C, sys.D);
Gwv = minreal(tf(num, den))

%% Bodeplot
h = figure(100)
bode(Gwv)
grid on
title('Transfer function from motor voltage to velocity')
saveas(h, 'motor to velocity.png');


%% Estimate transfer function for base system using LINEARIZE
% HOMEMADE from here down
% Velcity reference to tilt angle
load_system(model);
open_system(model);
% define points in model
ios(1) = linio(strcat(model,'/vel_ref'),1,'openinput');
ios(2) = linio(strcat(model, '/tilt_angle'),1,'openoutput');
% attach to model
setlinio(model,ios);
% Use the snapshot time(s) 0 seconds
op = [0];
% Linearize the model
sys = linearize(model,ios,op);
% get transfer function
[num,den] = ss2tf(sys.A, sys.B, sys.C, sys.D);
Gtv = minreal(tf(num, den))

%% Bodeplot
h = figure(100)
bode(Gtv)
grid on
title('Transfer function from reference velocity to tilt angle')
saveas(h, 'vel to tilt.png');

Kps = -1;              % From the Nyquist plot
t_ipost = 1 / 10;     % From bode plot of Gvt
save('control_variables.mat', 'Kps', 't_ipost', '-append');

%% Estimate transfer function for base system using LINEARIZE
% Transfer function Gtv_input to tilt_angle
% Velcity reference to tilt angle
load_system(model);
open_system(model);
% define points in model
ios(1) = linio(strcat(model,'/Gtv_input'),1,'openinput');
ios(2) = linio(strcat(model, '/tilt_angle'),1,'openoutput');
% attach to model
setlinio(model,ios);
% Use the snapshot time(s) 0 seconds
op = [0];
% Linearize the model
sys = linearize(model,ios,op);
% get transfer function
[num,den] = ss2tf(sys.A, sys.B, sys.C, sys.D);
Gtv_post = minreal(tf(num, den))

%% Bodeplot
h = figure(100)
bode(Gtv_post)
grid on
title('Transfer function from reference tilt angle to tilt angle')
saveas(h, 'tilt ref to tilt.png');

%% Tilt controller
Ni_tilt = 3;
alpha_tilt = 0.2;
phase_margin_tilt = 60;

% Find the new cross-over frequency
w_c_tilt = phaseBalance_equation(deg2rad(phase_margin_tilt), alpha_tilt, Ni_tilt, Gtv_post);

% Find the time constant ti and Cpi
ti_tilt = Ni_tilt / w_c_tilt;
Cpi_tilt = tf([ti_tilt, 1], [ti_tilt 0]);

% Find td and Cd
td_tilt = 1 / (w_c_tilt * sqrt(alpha_tilt));
Cd_tilt = tf([td_tilt, 1], [alpha_tilt*td_tilt, 1]);

% Find the Kp
open_loop_tf = Cpi_tilt * Cd_tilt * Gtv_post;
[mag, ~, ~] = bode(open_loop_tf, w_c_tilt); % 'bode' returns magnitude in absolute units (not dB), so no need to convert
Kp_tilt = 1 / squeeze(mag);
G_ol = Kp_tilt*Cpi_tilt*Cd_tilt*Gtv_post;

% save to matlab
save('control_variables.mat', 'Kp_tilt', '-append');

% This is a useful plot!
bode(G_ol);

% Calculate the closed loop tf when lead is in forward branch
G_cl_fwd = G_ol / (1 + G_ol);
% Calculate the closed loop tf when lead is in feedback branch
G_cl_fdb = (Kp_tilt*Cpi_tilt*Gtv_post) / (1 + Kp_tilt*Cpi_tilt*Gtv_post*Cd_tilt);

[num_Cd_tilt, den_Cd_tilt] = tfdata(Cd_tilt, 'v');
[num_Cpi_tilt, den_Cpi_tilt] = tfdata(Cpi_tilt, 'v');

% save to matlab
save('control_variables.mat', 'num_Cd_tilt', 'den_Cd_tilt', 'num_Cpi_tilt', 'den_Cpi_tilt', '-append');

figure;
step(G_cl_fwd);
hold on;
step(G_cl_fdb);
hold off;
legend('G_cl_fwd', 'G_cl_fdb');
disp(stepinfo(G_cl_fwd));
disp(stepinfo(G_cl_fdb));

%% VELOCITY CONTROLLER
% Transfer function from tilt reference to wheel velocity
% Velcity reference to tilt angle
load_system(model);
open_system(model);
% define points in model
ios(1) = linio(strcat(model,'/tilt_ref'),1,'openinput');
ios(2) = linio(strcat(model, '/wheel_vel_filter'),1,'openoutput');
% attach to model
setlinio(model,ios);
% Use the snapshot time(s) 0 seconds
op = [0];
% Linearize the model
sys = linearize(model,ios,op);
% get transfer function
[num,den] = ss2tf(sys.A, sys.B, sys.C, sys.D);
G_tilt_vel = minreal(tf(num, den)) % transfer function from tilt to velocity

%% Bodeplot
h = figure(100)
bode(G_tilt_vel)

grid on
title('Transfer function from tilt ref to wheel vel')
saveas(h, 'tilt ref to wheel vel.png');

%% Velcotiy controller parameters
% PI-lead-lag controller

Ni_vel = 5;
alpha_vel = 0.2;
phase_margin_vel = 60;
% since there was a chainsaw bode plot when just using PI-lead, the lag was
% also added
beta_vel = 10;

% calculate the PI lead controller
[wc_vel, Kp_vel, taui_vel, taud_vel, ok_vel] = findpidlag(G_tilt_vel, phase_margin_vel, Ni_vel, alpha_vel, beta_vel);
Cpi_vel = tf([taui_vel, 1], [taui_vel 1/beta_vel]);
Cd_vel = tf([taud_vel, 1], [alpha_vel*taud_vel, 1]);
G_ol_vel = Kp_vel*Cpi_vel*Cd_vel*G_tilt_vel;
G_cl_fwd_vel = G_ol_vel / (1 + G_ol_vel);
G_cl_fdb_vel = (Kp_vel*Cpi_vel*G_tilt_vel) / (1 + Kp_vel*Cpi_vel*G_tilt_vel*Cd_vel);

% get the individual numerator and denominators for the simulink model
[num_Cd_vel, den_Cd_vel] = tfdata(Cd_vel, 'v');
[num_Cpi_vel, den_Cpi_vel] = tfdata(Cpi_vel, 'v');

% save to matlab
save('control_variables.mat', 'num_Cd_vel', 'den_Cd_vel', 'num_Cpi_vel', 'den_Cpi_vel', '-append');
save('control_variables.mat', 'Kp_vel', '-append');

bode(G_ol_vel)

figure;
step(G_cl_fwd_vel);
hold on;
step(G_cl_fdb_vel);
hold off;
legend('G_cl_fwd', 'G_cl_fdb');
disp(stepinfo(G_cl_fwd_vel));
disp(stepinfo(G_cl_fdb_vel));

% Displaying the chainsaw effect.. 
% It is slowly oscilating around the given velocity (when velocity is zero)

%% Position controller
% Transfer function from velocity control reference to position
% Velcity controller reference to position
load_system(model);
open_system(model);
% define points in model
ios(1) = linio(strcat(model,'/vel_cont_ref'),1,'openinput');
ios(2) = linio(strcat(model, '/pos'),1,'openoutput');
% attach to model
setlinio(model,ios);
% Use the snapshot time(s) 0 seconds
op = [0];
% Linearize the model
sys = linearize(model,ios,op);
% get transfer function
[num,den] = ss2tf(sys.A, sys.B, sys.C, sys.D);
G_vel_pos = minreal(tf(num, den)) % transfer function from tilt to velocity

%% Bodeplot
h = figure(100)
bode(G_vel_pos)

grid on
title('Transfer function from velocity reference to position')
saveas(h, 'vel ref to position.png');

%% Position controller parameters
% will probably have to use a rate limiter or something similar as the
% signal towards velocity should be limited

Ni_pos = 5;
alpha_pos = 0.2;
phase_margin_pos = 60;

% calculate the PI lead controller
[wc_pos, Kp_pos, taui_pos, taud_pos, ok_pos] = findpid(G_vel_pos, phase_margin_pos, Ni_pos, alpha_pos);
Cpi_pos = tf([taui_pos, 1], [taui_pos 0]);
Cd_pos = tf([taud_pos, 1], [alpha_pos*taud_pos, 1]);
G_ol_pos = Kp_pos*Cpi_pos*Cd_pos*G_vel_pos;
G_cl_fwd_pos = G_ol_pos / (1 + G_ol_pos);
G_cl_fdb_pos = (Kp_pos*Cpi_pos*G_vel_pos) / (1 + Kp_pos*Cpi_pos*G_vel_pos*Cd_pos);

% get the individual numerator and denominators for the simulink model
[num_Cd_pos, den_Cd_pos] = tfdata(Cd_pos, 'v');
[num_Cpi_pos, den_Cpi_pos] = tfdata(Cpi_pos, 'v');

% save to matlab
save('control_variables.mat', 'num_Cd_pos', 'den_Cd_pos', 'num_Cpi_pos', 'den_Cpi_pos', '-append');
save('control_variables.mat', 'Kp_pos', '-append');

bode(G_ol_pos)

% the feedback position is more stable so going with that configuration
figure;
step(G_cl_fwd_pos);
hold on;
step(G_cl_fdb_pos);
hold off;
legend('G_cl_fwd', 'G_cl_fdb');
disp(stepinfo(G_cl_fwd_pos));
disp(stepinfo(G_cl_fdb_pos));