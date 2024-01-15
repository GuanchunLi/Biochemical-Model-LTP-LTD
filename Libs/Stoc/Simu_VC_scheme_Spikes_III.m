function [t_lst, y_VC_lst, n_channels_lst] = ...
    Simu_VC_scheme_Spikes_III(t_start, t_end, y_VC_0, n_channels_0, dt0, I_stim_func, NMDA_factor, CaV_factor, tk_pre_lst)

if length(tk_pre_lst) < 1
    error("Use the no-spike version simulation program!")
end

t_start_0 = t_start; t_end_0 = tk_pre_lst(1);
[t_lst, y_VC_lst, n_channels_lst] = ...
    Simu_VC_scheme_III(t_start_0, t_end_0, y_VC_0, n_channels_0, dt0, I_stim_func, NMDA_factor, CaV_factor);
n_channels = n_channels_lst(end, :)'; y_VC = y_VC_lst(end, :)';
% Adding the spikes
[n_channels] = Update_Channels_Spikes(n_channels);

for i = 1:(length(tk_pre_lst)-1)
    t_start_i = tk_pre_lst(i); t_end_i = tk_pre_lst(i+1);
    [t_lst_i, y_VC_lst_i, n_channels_lst_i] = ...
    Simu_VC_scheme_III(t_start_i, t_end_i, y_VC, n_channels, dt0, I_stim_func, NMDA_factor, CaV_factor);
    n_channels = n_channels_lst_i(end, :)'; y_VC = y_VC_lst_i(end, :)';
    t_lst = [t_lst; t_lst_i]; y_VC_lst = [y_VC_lst; y_VC_lst_i]; 
    n_channels_lst = [n_channels_lst; n_channels_lst_i];
    % Adding the spikes
    [n_channels] = Update_Channels_Spikes(n_channels);
end

t_start_i = tk_pre_lst(length(tk_pre_lst)); t_end_i = t_end;
[t_lst_i, y_VC_lst_i, n_channels_lst_i] = ...
Simu_VC_scheme_III(t_start_i, t_end_i, y_VC, n_channels, dt0, I_stim_func, NMDA_factor, CaV_factor);
t_lst = [t_lst; t_lst_i]; y_VC_lst = [y_VC_lst; y_VC_lst_i]; 
n_channels_lst = [n_channels_lst; n_channels_lst_i];

end