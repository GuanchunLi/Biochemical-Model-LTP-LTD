function [t_new, y_VC_new, n_channels_new, R_channels_new] = ...
    Simu_VC_nl_one_step_III(t, y_VC, n_channels, R_channels, I_stim, dt, NMDA_factor, CaV_factor)
% One step of stochastic simulation of nonlinear poisson process with
% time-varying rate (voltage & calcium dynamics)

% Compute matrix of dt (next reaction time)
% Compute the rates at t+dt
[y_VC_new_temp, dy_VC] = Update_VC_III(y_VC, n_channels, I_stim, dt, NMDA_factor, CaV_factor);
Vs_new = y_VC_new_temp(1);
c_Ca_SS_new = y_VC_new_temp(3);
[R_channels_new_temp] = Rate_Channels_III(t+dt, n_channels, Vs_new, c_Ca_SS_new);
% Generate the next reaction time
U = rand(size(R_channels));
Delta = R_channels.^2 - 2 * (R_channels_new_temp-R_channels) .* log(U) / dt;
T_channels = 0 * R_channels;
index = Delta > 0;
T_channels(~index) = 2*dt;
T_channels(index) = - 2 * log(U(index)) ./ (sqrt(Delta(index)) + R_channels(index));

% Update the system
[t_new, y_VC_new, n_channels_new, R_channels_new] = ...
Update_VC_family_III(t, y_VC, n_channels, dt, ...
R_channels_new_temp, T_channels, y_VC_new_temp, dy_VC);
end