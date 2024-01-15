function [t_lst, y_lst] = Simu_stoch_scheme(T, y0, Reaction_stoc, Update_stoc, paras, cal_source)
% SIMU_STOCH_SCHEME
% The scheme used to simulate the CaMKII reactions (Stochastic version)
t = 0; y = y0;
k = 1;
t_lst(k,:) = t; y_lst(k, :, :) = y;
while t < T
    [t_new, idx] = Reaction_stoc(t, y, cal_source(t), paras.c_E0, paras.c_CaM, paras);
    y = Update_stoc(y, idx, paras);
    t = t_new;
    % fprintf(['t = ', num2str(t), ' n_MP_f = ', num2str(y(8)), ' n_C = ', num2str(y(9)), '\n']);
    k = k+1;
    t_lst(k) = t;
    y_lst(k, :, :) = y;
    if mod(k, 100) == 0
        fprintf(sprintf('t = %f, step = %d \n', t, k));
    end
end