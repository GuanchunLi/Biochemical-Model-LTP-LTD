function [state_final] = Simu_Stoc_scheme_norecord(t_start, t_end, dt0, state_0, paras, cal_source)

% Model Set Up
Rate_func = @(t, y) Rate_Stoc_CaMKII_family(t,y, cal_source(t), paras);
% Rate_func = @(t, y) Rate_Stoc_CaMKII_family_masked(t,y, cal_source(t), paras);
Update_func = @(t, y, Rate_new_temp, T, dt, Rate_func) Update_Stoc_CaMKII_family(t, y, Rate_new_temp, T, dt, Rate_func, paras);

% Initilize
t = t_start; state = state_0; k = 1;
% N_timestep = 200000;
% state_lst = zeros([size(state), N_timestep]); state_lst(:, :, k) = state;
% t_lst = zeros(N_timestep, 1); t_lst(k) = t;
[R_mat_all] = Rate_Stoc_CaMKII_family(t, state, cal_source(t), paras);
% Rate_lst = zeros([size(R_mat_all), N_timestep]); Rate_lst(:, :, k) = R_mat_all;

% Run the simulation
while t < t_end
    dt = min(t_end-t, dt0);
    [t_new, state_new, Rate_new] = Simu_stoch_nl_one_step(t, state, R_mat_all, ...
        Rate_func, Update_func, dt);
    t = t_new; state = state_new; R_mat_all = Rate_new;
%     k = k + 1; state_lst(:, :, k) = state; t_lst(k) = t;
    % Rate_lst(:, :, k) = R_mat_all;
end

state_final = state;
% t_lst = t_lst(1:k);
% state_lst = state_lst(:, :, 1:k);
% Rate_lst = Rate_lst(:, :, 1:k);

end