function [t_new, y_VC_new, n_channels_new, R_channels_new] = ...
Update_VC_family(t, y_VC, n_channels, dt, ...
R_channels_new_temp, T_channels, y_VC_new_temp, dy_VC)

[dt_update, idx_update] = min(T_channels);

if dt_update > dt % No reaction occurs
    t_new = t + dt;
    y_VC_new = y_VC_new_temp;
    n_channels_new = n_channels;
    R_channels_new = R_channels_new_temp;

elseif dt_update > 0 % Reaction occurs
    t_new = t + dt_update;
    [n_channels_new] = Update_Channels(n_channels, idx_update);
    [y_VC_new] = y_VC + dy_VC * dt_update;
    Vs_new = y_VC_new(1);
    c_Ca_SS_new = y_VC_new(3);
    [R_channels_new] = Rate_Channels(t_new, n_channels_new, Vs_new, c_Ca_SS_new);

else
    error("dt < 0!");
end

end