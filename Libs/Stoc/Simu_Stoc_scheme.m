function [t_lst, state, c_Ms_lst, c_MPs_lst, c_E0_lst] = Simu_Stoc_scheme(t_start, t_end, dt0, state_0, paras, cal_source)

state_save_flag = 0; % paras.state_save_flag; % 0; % no state_lst

% Model Set Up
Rate_func = @(t, y) Rate_Stoc_CaMKII_family(t,y, cal_source(t), paras);
% Rate_func = @(t, y) Rate_Stoc_CaMKII_family_masked(t,y, cal_source(t), paras);
Update_func = @(t, y, Rate_new_temp, T, dt, Rate_func) Update_Stoc_CaMKII_family(t, y, Rate_new_temp, T, dt, Rate_func, paras);

% Initilize
t = t_start; state = state_0; k = 1;
[c_Ms, c_MPs] = State2Con(state, paras); c_E0 = state(1, end) / paras.coef_V;
N_timestep = 100000;
% if state_save_flag
%     state_lst = zeros([size(state), N_timestep]); state_lst(:, :, k) = state;
% end
t_lst = zeros(N_timestep, 1); t_lst(k) = t;
c_Ms_lst = zeros(N_timestep, 1); c_Ms_lst(k) = c_Ms;
c_MPs_lst = zeros(N_timestep, 1); c_MPs_lst(k) = c_MPs;
c_E0_lst = zeros(N_timestep, 1); c_E0_lst(k) = c_E0;
[R_mat_all] = Rate_Stoc_CaMKII_family(t, state, cal_source(t), paras);
% Rate_lst = zeros([size(R_mat_all), N_timestep]); Rate_lst(:, :, k) = R_mat_all;

% Run the simulation
while t < t_end
    dt = min(t_end-t, dt0);
    [t_new, state_new, Rate_new] = Simu_stoch_nl_one_step(t, state, R_mat_all, ...
        Rate_func, Update_func, dt);
    t = t_new; state = state_new; R_mat_all = Rate_new;
    k = k + 1; t_lst(k) = t; 
%     if state_save_flag
%        state_lst(:, :, k) = state; 
%     end
    [c_Ms, c_MPs] = State2Con(state, paras); c_E0 = state(1, end) / paras.coef_V;
    c_Ms_lst(k) = c_Ms; c_MPs_lst(k) = c_MPs; c_E0_lst(k) = c_E0;
    % Rate_lst(:, :, k) = R_mat_all;
end

t_lst = t_lst(1:k);
% state_lst = state_lst(:, :, 1:k);
% Rate_lst = Rate_lst(:, :, 1:k);
c_Ms_lst = c_Ms_lst(1:k);
c_MPs_lst = c_MPs_lst(1:k);
c_E0_lst = c_E0_lst(1:k);

end