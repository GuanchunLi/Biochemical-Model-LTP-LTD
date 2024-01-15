function [t_lst, c_Ms_lst, c_MPs_lst, c_E0_lst, state] = Simu_Stoc_scheme_multistep(t_start, t_end, ...
          t_cri_lst, max_step, d_step, state_0, paras, cal_source)

state_save_flag = 0; % paras.state_save_flag; % 0 % no state_lst

t_cri_lst = t_cri_lst((t_cri_lst >= t_start) .* (t_cri_lst <= t_end) > 0);
n_cri = length(t_cri_lst);
t_simu_lst = zeros(2*n_cri, 1);
max_step_lst = zeros(2*n_cri, 1);

if isempty(t_cri_lst)
    [t_lst, state, c_Ms_lst, c_MPs_lst, c_E0_lst] = ...
    Simu_Stoc_scheme(t_start, t_end, 1e-3, state_0, paras, cal_source);
else
% Initialize
t_simu_lst(1) = t_cri_lst(1) - d_step;
max_step_lst(1) = 1;

k = 1; i = 1;
dt_lst = diff(t_cri_lst);

while i <= (n_cri-1)
    if dt_lst(i) > 2*d_step
        k = k + 1;
        tk_end = t_cri_lst(i) + d_step;
        t_simu_lst(k) = tk_end; max_step_lst(k) = 1e-4;
        k = k + 1;
        tk_end = t_cri_lst(i+1) - d_step;
        t_simu_lst(k) = tk_end; max_step_lst(k) = max_step;
        i = i + 1;
    else
        k = k + 1;
        tk_end = t_cri_lst(i) + d_step / 2;
        t_simu_lst(k) = tk_end; max_step_lst(k) = 1e-4;
        i = i + 1;
    end
end

k = k + 1;
tk_end = t_cri_lst(n_cri) + d_step; 
t_simu_lst(k) = tk_end; max_step_lst(k) = 1e-4;

k = k + 1;
t_simu_lst(k) = t_end; max_step_lst(k) = 1;

t_simu_lst =  t_simu_lst(1:k); max_step_lst = max_step_lst(1:k);

% if state_save_flag
%     state_lst = [];
% end

% Start the simulation
[t_lst, state, c_Ms_lst, c_MPs_lst, c_E0_lst] = Simu_Stoc_scheme(t_start, t_simu_lst(1), ...
     max_step_lst(1), state_0, paras, cal_source);
% state_lst = cat(3, state_lst, state); state = squeeze(state(:, :, end));

for i = 1:(length(t_simu_lst)-1)
    [t_lst_i, state, c_Ms_lst_i, c_MPs_lst_i, c_E0_lst_i] = Simu_Stoc_scheme(t_simu_lst(i), t_simu_lst(i+1), ...
        max_step_lst(i+1), state, paras, cal_source);
    t_lst = [t_lst; t_lst_i]; c_E0_lst = [c_E0_lst; c_E0_lst_i];
    c_Ms_lst = [c_Ms_lst; c_Ms_lst_i]; c_MPs_lst = [c_MPs_lst; c_MPs_lst_i];
%     if state_save_flag
%         state_lst = cat(3, state_lst, state); state = squeeze(state(:, :, end));
%     end
end

end
