function [t_lst, y_VC_lst, n_channels_lst] = ...
    Simu_VC_scheme_III(t_start, t_end, y_VC_0, n_channels_0, dt0, I_stim_func, NMDA_factor, CaV_factor)
% SIMU_STOCH_SCHEME
% The scheme used to simulate the CaMKII reactions (Stochastic version)

% Initilize
t = t_start; y_VC = y_VC_0; n_channels = n_channels_0; k = 1;
% How large for this term?
N_timestep = 100000 * 180;
t_lst = zeros(N_timestep, 1); t_lst(k) = t;
y_VC_lst = zeros([N_timestep, length(y_VC)]); y_VC_lst(k, :) = y_VC';
n_channels_lst = zeros([N_timestep, length(n_channels)]); n_channels_lst(k, :) = n_channels';
[R_channels] = Rate_Channels_III(t, n_channels, y_VC(1), y_VC(3));

while t < t_end
    dt = min(t_end-t, dt0);
    [t_new, y_VC_new, n_channels_new, R_channels_new] = ...
    Simu_VC_nl_one_step_III(t, y_VC, n_channels, R_channels, I_stim_func(t), dt, NMDA_factor(t), CaV_factor(t));
    t = t_new; y_VC = y_VC_new; n_channels = n_channels_new; R_channels = R_channels_new;
    k = k + 1; t_lst(k) = t; y_VC_lst(k, :) = y_VC'; n_channels_lst(k, :) = n_channels';
end

t_lst = t_lst(1:k);
y_VC_lst = y_VC_lst(1:k, :);
n_channels_lst = n_channels_lst(1:k, :);

end
