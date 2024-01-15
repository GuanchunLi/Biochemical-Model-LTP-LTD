function [t_Ca_lst, Ca_lst, t_lst, c_MPs_lst, c_E0_lst, y_lst] = ...
Simu_Framework_III(t_end, tk_pre_lst, I_stim_func, t_cri_lst, paras, ...
    stoc_input_flag, stoc_reaction_flag, save_state_flag)
%UNTITLED 此处提供此函数的摘要

% y_lst -- for determinstic, all CaMKII states; for stochastic, state_lst
%% Simulate Calcium input
if stoc_input_flag % Stochastic calcium input
    % Initialize
    NMDA_factor = @(t) paras.NMDA_factor * 1.305;
    CaV_factor = @(t) paras.CaV_factor;
    dt0 = 1e-5;
    y_VC_0 = [-71.4479454884918; 0.103233668741414; 0.122470478619358; 0.00756800588549419; 0.980292454270282; 0.160041702948184];
    y_VC_0(7:12) = paras.y0_Ca(4:9);  % Variables of CaV12 model
    y_VC_0(13:18) = paras.y0_Ca(10:15); % Variables of CaV13 model
    n_channels_0 = Con2State_VC_III(paras.y0_Ca, paras);

    [t_lst, y_VC_lst, ~] = ...
        Simu_VC_scheme_Spikes_III(0, t_end, y_VC_0, n_channels_0, dt0, I_stim_func, NMDA_factor, CaV_factor, tk_pre_lst);
    t_Ca_lst = t_lst; Ca_lst = y_VC_lst(:, 2);
    [t_Ca_lst, t_index] = unique(t_Ca_lst);
    Ca_lst = Ca_lst(t_index);
else % Determinstic calcium input
    % Initialize
    NMDA_factor = @(t) paras.NMDA_factor;
    CaV_factor = @(t) paras.CaV_factor;
    tau_x_AMPA = 0.05/(10^3); alpha_x_AMPA = 1; 
    x_AMPA_func = @(t) x_spike(t, tk_pre_lst, tau_x_AMPA, alpha_x_AMPA);
    tau_x_NMDA = 2/(10^3); alpha_x_NMDA = 1;
    x_NMDA_func = @(t) x_spike(t, tk_pre_lst, tau_x_NMDA, alpha_x_NMDA);

    opts = odeset('RelTol',1e-7, 'AbsTol',1e-6, 'MaxStep',1e-3);
    circuit_model = @(t, y) V_dynamic_Brunel_stoc(y, x_AMPA_func(t),...
                               x_NMDA_func(t), I_stim_func(t), NMDA_factor(t), CaV_factor(t));
    [t_dc, y_dc] = ode15s(@(t, y) V_dynamic_Brunel_stoc(y, 0, 0, 0, NMDA_factor(t), CaV_factor(t)), [0,10], paras.y0_Ca, opts);
    % y0 = paras.y0_Ca;
    y0 = y_dc(end,:)';
    sol_Ca = ode15s(circuit_model, [0, t_end], y0, opts);
    t = sol_Ca.x; y = sol_Ca.y';
    t_Ca_lst = t; Ca_lst = y(:, 16);
end

%% Simulate CaMKII reactions
if stoc_reaction_flag % Stochastic reaction dynamics
    % Initialize
    max_step = 1e-3; d_step = 1e-3;
    state = Con2State_CaMKII_family(paras, paras.y0_CaMKII_simu);
    
    % Start the simulation
    T = t_end; n_range = 32;
    T_lst = linspace(0, T, n_range);
    t_lst = []; y_lst = []; 
    c_E0_lst = []; c_Ms_lst = []; c_MPs_lst = [];
    for i = 1:(n_range-1)
        t0 = T_lst(i); t1 = T_lst(i+1);
        t_index = (t_Ca_lst >= (t0-1e-2)) .* (t_Ca_lst <= (t1+1e-2)) > 0;
        t_Ca_lst_i = t_Ca_lst(t_index);
        Ca_lst_i = Ca_lst(t_index);
        cal_source_i = @(t) interp1(t_Ca_lst_i, Ca_lst_i, t);
        [t_lst_i, c_Ms_lst_i, c_MPs_lst_i, c_E0_lst_i, state] = Simu_Stoc_scheme_multistep(t0, t1, t_cri_lst, max_step, d_step, state, paras, cal_source_i);
        t_lst = [t_lst; t_lst_i]; c_E0_lst = [c_E0_lst; c_E0_lst_i];
        c_Ms_lst = [c_Ms_lst; c_Ms_lst_i]; c_MPs_lst = [c_MPs_lst; c_MPs_lst_i];
    end
    % t_cri_lst = unique(sort([tk_pre_lst, tk_post_lst]));
    % [t_lst, state_lst] = Simu_Stoc_scheme_multistep(0, T, t_cri_lst, max_step, d_step, state_0, paras, cal_source);    
    y_lst = [];

else % Determinstic reaction dynamics
    % t_cri_lst = unique(sort([tk_pre_lst, tk_post_lst]));
    y = paras.y0_CaMKII_simu;
    T = t_end; n_range = 32;
    T_lst = linspace(0, T, n_range);
    t_lst = []; y_lst = [];
    for i = 1:(n_range-1)
        t0 = T_lst(i); t1 = T_lst(i+1);
        t_index = (t_Ca_lst >= (t0-1e-2)) .* (t_Ca_lst <= (t1+1e-2)) > 0;
        t_Ca_lst_i = t_Ca_lst(t_index);
        Ca_lst_i = Ca_lst(t_index);
        cal_source_i = @(t) interp1(t_Ca_lst_i, Ca_lst_i, t);
        [t_lst_i, y_lst_i] = Simu_ode_scheme_multistep(t0, t1, t_cri_lst, 1, 1e-3, y, paras, cal_source_i);
        t_lst = [t_lst, t_lst_i]; y_lst = [y_lst; y_lst_i]; y = y_lst_i(end, :)';
    end
    c_MPs_lst = y_lst(:, 6:19) * paras.MP_lst';
    c_E0_lst= y_lst(:, 1);
end

end