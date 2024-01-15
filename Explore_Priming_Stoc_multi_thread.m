%% Initialization
clear;

addpath(genpath([pwd, '/Libs/']));
addpath(genpath([pwd, '/Parameters']));

% Set up all parameters for the model
paras.n = 6;    % Number of monomers in each CaMKII (consider only one ring)
paras.v = 1;    % Reaction speed scale (reaction_rate = base_rate * v)
% paras = Param_CaM(paras);
paras.v_tot = 2;
paras = Param_CaM(paras);
paras = Param_Concent_Retuning(paras);
paras = Param_FastEqual(paras);
paras = Param_ReactRate_Retuning(paras);
paras = Param_PhosTrans(paras);
paras = Param_PP1_mod3_Retuning(paras);
paras = Param_RateModel(paras);
paras = Param_Readout(paras);
paras = Param_InitVal_Retuning(paras);
paras.display = 0;

%% Generate Calcium source
t_end = 120;
tk_pre_lst = [15:0.1:15.2]; % Pre-spike time list, unit: s

parpool(32);

flag = 6;
% dc_lst = 0:1:350;
dc_lst = linspace(0, 384, 256);
N_dc = length(dc_lst);
N_trail = 100;
CaMKII_mat = zeros(N_trail, N_dc);
E_mat = zeros(N_trail, N_dc);
AMPA_mat = zeros(N_trail, N_dc);
V_mat = zeros(N_trail, N_dc);

parfor j = 1:length(dc_lst)
    dc = dc_lst(j);
    for k = 1:N_trail
        [CaMKII_SS, E_SS, AMPA_SS, E_base] = run_simulation_k(tk_pre_lst, t_end, dc, k, paras, flag);
        CaMKII_mat(k, j) = CaMKII_SS;
        E_mat(k, j) = E_SS;
        AMPA_mat(k, j) = AMPA_SS;
        V_mat(k, j) = E_base;
    end
    fprintf('%d finished! \n', j);
end

fig1 = figure(1);
plot(dc_lst, mean(CaMKII_mat));
hold on;
plot(dc_lst, mean(E_mat));
saveas(fig1, sprintf('Priming_Dependency_Concent_flag=%d.png', flag));
% 
fig2 = figure(2);
plot(dc_lst, mean(AMPA_mat-1));
saveas(fig2, sprintf('Priming_Dependency_AMPA_flag=%d.png', flag));

save(sprintf('Priming_Dependency_flag=%d.mat', flag), 'dc_lst', 'CaMKII_mat', 'E_mat', 'AMPA_mat', 'V_mat');

function [CaMKII_SS, E_SS, AMPA_SS, E_base] = run_simulation_k(tk_pre_lst, t_end, dc, k, paras, flag)
sigma_NMDA = 3.3/100;
NMDA_factor_lst = 0 * tk_pre_lst + 1;
NMDA_factor_lst = max(NMDA_factor_lst + sigma_NMDA * randn(1, length(tk_pre_lst)),0);
NMDA_factor = @(t) interp1([0,tk_pre_lst, t_end], [1,NMDA_factor_lst,1], t, 'previous');

tk_post_lst = max(tk_pre_lst) + 500/1000;
sigma_CaV = 10/100;
CaV_factor_lst = 0 * tk_post_lst + 1;
CaV_factor_lst = max(CaV_factor_lst + sigma_CaV * randn(1, length(tk_post_lst)),0);
CaV_factor = @(t) interp1([0,tk_post_lst, max(tk_post_lst)+1, t_end], [1,CaV_factor_lst,1,1], t, 'previous');

tau_x_AMPA = 0.05/(10^3); alpha_x_AMPA = 1; 
x_AMPA_func = @(t) x_spike(t, tk_pre_lst, tau_x_AMPA, alpha_x_AMPA);
tau_x_NMDA = 2/(10^3); alpha_x_NMDA = 1;
x_NMDA_func = @(t) x_spike(t, tk_pre_lst, tau_x_NMDA, alpha_x_NMDA);
duration = 300 + 2;
if flag == 1 % Control case
    I_stim_func = @(t) 3000 * any((tk_post_lst < t-2/1e3) .* (tk_post_lst > t-duration/1e3)) + dc; % * (t>tk_post_lst);
elseif flag == 2 % Depression of CaV1.3 for all times
    I_stim_func = @(t) 3000 * any((tk_post_lst < t-2/1e3) .* (tk_post_lst > t-duration/1e3)) + dc * 0;
elseif flag == 3 % Depression of CaV1.3 except plateau potential
    I_stim_func = @(t) (3000+dc) * any((tk_post_lst < t-2/1e3) .* (tk_post_lst > t-duration/1e3));
elseif flag == 4 % Depression of CaV1.3 only during ‘Priming’ (before plateau potential)
    I_stim_func = @(t) 3000 * any((tk_post_lst < t-2/1e3) .* (tk_post_lst > t-duration/1e3)) + dc * (t>tk_post_lst);
elseif flag == 5 % Depression of CaV1.3 only during ‘Priming’ (before plateau potential)
    I_stim_func = @(t) 3000 * any((tk_post_lst < t-2/1e3) .* (tk_post_lst > t-duration/1e3)) + dc * (t<min(tk_pre_lst));
elseif flag == 6 % Depression of CaV1.3 only except pre-spikes
    I_stim_func = @(t) 3000 * any((tk_post_lst < t-2/1e3) .* (tk_post_lst > t-duration/1e3)) ...
               + dc * (t>min(tk_pre_lst)) * (t<max(tk_pre_lst));
end

circuit_model = @(t, y) V_dynamic_Brunel_stoc(y, x_AMPA_func(t),...
                           x_NMDA_func(t), I_stim_func(t),...
                           ... % 0.412, 2.02);
                           0.412*NMDA_factor(t), 2.02*CaV_factor(t));

opts = odeset('RelTol',1e-7, 'AbsTol',1e-6, 'MaxStep',1e-3);
[t_dc, y_dc] = ode15s(@(t, y) V_dynamic_Brunel_stoc(y, 0, 0, I_stim_func(t), 0.412, 2.02), [0,3], paras.y0_Ca, opts);
y0 = y_dc(end,:)';

% y0 = paras.y0_Ca;

sol_Ca = ode15s(circuit_model, [0, t_end], y0, opts);
t = sol_Ca.x; y = sol_Ca.y';
E_base = y(end, 1);
cal_source = @(t) Ca_func(sol_Ca, t);

%% Initial state of CaMKII system

% Baseline SS
y0 = paras.y0_CaMKII;

%% Start the simulation
y = y0; T = t_end;
t_cri_lst = unique(sort([tk_pre_lst, tk_post_lst]));
[t_lst, y_lst] = Simu_ode_scheme_multistep(0, T, t_cri_lst, 1, 1e-3, y, paras, cal_source);

y_ss = y_lst(end, :);
CaMKII_SS = y_ss(6:19) * paras.MP_lst' / (paras.n * paras.c0) * 100;
E_SS = y_ss(1); 
AMPA_SS = y_ss(20) / y0(20);

% fprintf('%d finished! \n', k)
end

function [c_Ca] = Ca_func(sol, t)
y = deval(sol, t);
c_Ca = y(16);
end
