function [t_lst, y_lst] = Simu_ode_scheme_multistep(t_start, t_end, ...
          t_cri_lst, max_step, d_step, y0, paras, cal_source)

t_cri_lst = t_cri_lst((t_cri_lst >= t_start) .* (t_cri_lst <= t_end) > 0);
n_cri = length(t_cri_lst);
t_simu_lst = zeros(2*n_cri, 1);
max_step_lst = zeros(2*n_cri, 1);

if isempty(t_cri_lst)
    [t_lst, y_lst] = Simu_ode_scheme(t_start, t_end, ...
    min(1e-1, (t_end-t_start)/10), y0, paras, cal_source);
else
% Initialize
t_simu_lst(1) = t_cri_lst(1) - d_step;
max_step_lst(1) = max_step;

k = 1; i = 1;
dt_lst = diff(t_cri_lst);

while i <= (n_cri-1)
    if dt_lst(i) > 2*d_step
        k = k + 1;
        tk_end = t_cri_lst(i) + d_step;
        t_simu_lst(k) = tk_end; max_step_lst(k) = 1e-1;
        k = k + 1;
        tk_end = t_cri_lst(i+1) - d_step;
        t_simu_lst(k) = tk_end; max_step_lst(k) = max_step;
        i = i + 1;
    else
        k = k + 1;
        tk_end = t_cri_lst(i) + d_step / 2;
        t_simu_lst(k) = tk_end; max_step_lst(k) = 1e-1;
        i = i + 1;
    end
end

k = k + 1;
tk_end = t_cri_lst(n_cri) + d_step; 
t_simu_lst(k) = tk_end; max_step_lst(k) = 1e-1;

k = k + 1;
t_simu_lst(k) = t_end; max_step_lst(k) = max_step;

t_simu_lst =  t_simu_lst(1:k); max_step_lst = max_step_lst(1:k);
[t_lst, y_lst] = Simu_ode_scheme(t_start, t_simu_lst(1), ...
     max_step_lst(1), y0, paras, cal_source);
y = y_lst(end, :)';
for i = 1:(length(t_simu_lst)-1)
    [t_lst_i, y_lst_i] = Simu_ode_scheme(t_simu_lst(i), t_simu_lst(i+1), ...
        max_step_lst(i+1), y, paras, cal_source);
    t_lst = [t_lst, t_lst_i]; y_lst = [y_lst; y_lst_i]; y = y_lst_i(end, :)';
end
end
end
% % Initilization
% [t_lst, y_lst] = Simu_ode_scheme(t_start, t_cri_lst(1) - d_step, ...
%     max_step_lst(1), y0, paras, cal_source);
% y = y_lst(end, :)';
% 
% N = length(t_cri_lst);
% for i = 1:(N)
%     t_start_i = max(t_start, t_cri_lst(i) - d_step);
%     t_end_i = min(t_end, t_cri_lst(i) + d_step);
%     [t_lst_i, y_lst_i] = Simu_ode_scheme(t_start_i, t_end_i, max_step_lst(2*i), ...
%                                      y, paras, cal_source);
%     t_lst = [t_lst, t_lst_i]; y_lst = [y_lst; y_lst_i]; y = y_lst_i(end, :)';
%     t_start_j = t_end_i;
%     if i < N
%         t_end_j = min(t_end, t_cri_lst(i+1) - d_step);
%     else
%         t_end_j = t_end;
%     end
%     [t_lst_j, y_lst_j] = Simu_ode_scheme(t_start_j, t_end_j, max_step_lst(2*i+1), ...
%                                      y, paras, cal_source);
%     t_lst = [t_lst, t_lst_j]; y_lst = [y_lst; y_lst_j]; y = y_lst_j(end, :)';
% end