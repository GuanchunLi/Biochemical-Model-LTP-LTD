function [state_final] = Simu_Stoc_scheme_multistep_norecord(t_start, t_end, ...
          t_cri_lst, max_step, d_step, state_0, paras, cal_source)

n_cri = length(t_cri_lst);
t_simu_lst = zeros(2*n_cri, 1);
max_step_lst = zeros(2*n_cri, 1);

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

% Start the simulation
[state] = Simu_Stoc_scheme_norecord(t_start, t_simu_lst(1), ...
     max_step_lst(1), state_0, paras, cal_source);
% state = squeeze(state_lst(:, :, end));
for i = 1:(length(t_simu_lst)-1)
    [state] = Simu_Stoc_scheme_norecord(t_simu_lst(i), t_simu_lst(i+1), ...
        max_step_lst(i+1), state, paras, cal_source);
%     t_lst = [t_lst; t_lst_i]; state_lst = cat(3, state_lst, state_lst_i); state = squeeze(state_lst(:, :, end));
end
state_final = state;
end
