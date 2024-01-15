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
paras.c0 = sum(paras.y0_CaMKII(5:19));

%% Priming Dynamics for different scenarios
t_end = 120;
tk_pre_lst = [15:0.1:15.2]; % Pre-spike time list, unit: s

flag_lst = [7];
dc = 90; flag_CaV = 1;

for k = 1:length(flag_lst)
    flag = flag_lst(k);
    [t_Ca_lst, Ca_lst, V_lst, t_lst, y_lst, c_MP_lst, c_M_lst, c_E_lst] = run_simulation_k(tk_pre_lst, t_end, dc, k, paras, flag, flag_CaV);
    save(sprintf('Result/Result_04_28/Priming/Priming_dynamic_dc90_flag=%d.mat', flag), 'flag', 'dc', ...
     't_Ca_lst', 'Ca_lst', 't_lst', 'y_lst', 'c_MP_lst', 'c_E_lst', 'c_M_lst');
    fig = figure();
    plot(t_lst, 100*c_MP_lst/(paras.n*paras.c0));
    hold on;
    plot(t_lst, 100*c_M_lst/(paras.n*paras.c0));
    plot(t_lst, 100*c_E_lst/paras.c_E_tot)
    xlabel('Time (s) ', 'FontSize', 18);
    ylabel('Activity (%)', 'FontSize', 18);
    legend('pCaMKII', 'oCaMKII', 'PP1');
    saveas(fig, sprintf('Result/Result_04_28/Priming/Priming_dynamic_dc90_flag=%d.png', flag));
end

%% Priming Dynamics for different DC currents
% t_end = 120;
% tk_pre_lst = [15:0.1:15.2]; % Pre-spike time list, unit: s
% % tk_pre_lst = [];
% % 
% % flag = [4,5,8];
% flag = 1; flag_CaV = 1;
% dc_lst = [90];
% 
% for k = 1:length(dc_lst)
%     dc = dc_lst(k);
%     [t_Ca_lst, Ca_lst, V_lst, t_lst, y_lst, c_MP_lst, c_M_lst, c_E_lst] = run_simulation_k(tk_pre_lst, t_end, dc, k, paras, flag, flag_CaV);
%     save('Result/Figures_For_Paper/Fig4/Priming_dynamic_dc90_NMDAR_Block.mat', 'flag', 'dc', ...
%      't_Ca_lst', 'Ca_lst', 't_lst', 'y_lst', 'c_MP_lst', 'c_E_lst', 'c_M_lst');
%     fig = figure();
%     plot(t_lst, 100*c_MP_lst/(paras.n*paras.c0));
%     hold on;
%     plot(t_lst, 100*c_M_lst/(paras.n*paras.c0));
%     plot(t_lst, 100*c_E_lst/paras.c_E_tot)
%     xlabel('Time (s) ', 'FontSize', 18);
%     ylabel('Activity (%)', 'FontSize', 18);
%     legend('pCaMKII', 'oCaMKII', 'PP1');
%     saveas(fig, 'Result/Figures_For_Paper/Fig4/Priming_dynamic_dc90_NMDAR_Block.png');
% end       

%% Aux functions
function [t_Ca_lst, Ca_lst, V_lst, t_lst, y_lst, c_MP_lst, c_M_lst, c_E_lst] = run_simulation_k(tk_pre_lst, t_end, dc, k, paras, flag, flag_CaV)
sigma_NMDA = 3.3/100;
NMDA_factor_lst = 0 * tk_pre_lst + 1;
% NMDA_factor_lst = max(NMDA_factor_lst + sigma_NMDA * randn(1, length(tk_pre_lst)),0);
NMDA_factor = @(t) interp1([0,tk_pre_lst, t_end], [1,NMDA_factor_lst,1], t, 'previous');

if isempty(tk_pre_lst)
    tk_post_lst = 40;
else
    tk_post_lst = max(tk_pre_lst) + 500/1000;
end

sigma_CaV = 10/100;
CaV_factor_lst = 0 * tk_post_lst + 1;
% CaV_factor_lst = max(CaV_factor_lst + sigma_CaV * randn(1, length(tk_post_lst)),0);
CaV_factor = @(t) interp1([0,tk_post_lst, max(tk_post_lst)+1e-1, t_end], [1,CaV_factor_lst,1,1], t, 'previous');

% CaV_F = 2.02;
% CaV_factor = @(t) CaV_F;
tau_x_AMPA = 0.05/(10^3); alpha_x_AMPA = 1; 
x_AMPA_func = @(t) x_spike(t, tk_pre_lst, tau_x_AMPA, alpha_x_AMPA);
tau_x_NMDA = 2/(10^3); alpha_x_NMDA = 1;
x_NMDA_func = @(t) x_spike(t, tk_pre_lst, tau_x_NMDA, alpha_x_NMDA);
duration = 300 + 2;
if flag == 1 % Control case
    I_stim_func = @(t) 4000 * any((tk_post_lst < t-2/1e3) .* (tk_post_lst > t-duration/1e3)) + dc; % * (t>tk_post_lst);
elseif flag == 2 % Depression of CaV1.3 for all times
    I_stim_func = @(t) 4000 * any((tk_post_lst < t-2/1e3) .* (tk_post_lst > t-duration/1e3)) + dc * 0;
elseif flag == 3 % Depression of CaV1.3 except plateau potential
    I_stim_func = @(t) (4000+dc) * any((tk_post_lst < t-2/1e3) .* (tk_post_lst > t-duration/1e3));
elseif flag == 4 % DC after plateau potential
    I_stim_func = @(t) 4000 * any((tk_post_lst < t-2/1e3) .* (tk_post_lst > t-duration/1e3)) + dc * (t>tk_post_lst);
elseif flag == 5 % DC before pre-synaptic spikes
    I_stim_func = @(t) 4000 * any((tk_post_lst < t-2/1e3) .* (tk_post_lst > t-duration/1e3)) + dc * (t<min(tk_pre_lst));
elseif flag == 7 % DC only during pre-synaptic spikes
    I_stim_func = @(t) 3000 * any((tk_post_lst < t-2/1e3) .* (tk_post_lst > t-duration/1e3)) ...
        + dc * ((t>min(tk_pre_lst)) .* (t<max(tk_pre_lst)));
elseif flag == 8 % DC before plateau potential
    I_stim_func = @(t) 4000 * any((tk_post_lst < t-2/1e3) .* (tk_post_lst > t-duration/1e3)) ...
        + dc * (t<tk_post_lst);
end

CaV12_factor = CaV_factor;
CaV13_factor = CaV_factor;

if flag_CaV == 0 % Shut down CaV1.3 during priming
    CaV13_factor = @(t) CaV_factor(t) * (t>min(tk_post_lst));
end

opts = odeset('RelTol',1e-7, 'AbsTol',1e-6, 'MaxStep',1e-3);
circuit_model = @(t, y) V_dynamic_Brunel_stoc(y, x_AMPA_func(t),...
                           x_NMDA_func(t), I_stim_func(t), ... % 0.412, 2.02); % ...
                           0.412*NMDA_factor(t), 2.02*CaV_factor(t));
dc_model = @(t, y) V_dynamic_Brunel_stoc(y, 0, 0, dc, 0.412*NMDA_factor(t), 2.02*CaV_factor(t));
[t_dc, y_dc] = ode15s(dc_model, [0, 5], paras.y0_Ca, opts);
% y0 = paras.y0_Ca;
y0 = y_dc(end,:)';

sol_Ca = ode15s(circuit_model, [0, t_end], y0, opts);
t = sol_Ca.x; y = sol_Ca.y';
t_Ca_lst = t; Ca_lst = y(:, 16); V_lst = y(:, 1);
% E_base = y(end, 1);
cal_source = @(t) Ca_func(sol_Ca, t);

%% Initial state of CaMKII system

% Baseline SS
y0 = paras.y0_CaMKII;
% % y0(5:19) = y0(5:19) * 1.5; 
% y0(5) = y0(5) + 6;
paras.c0 = sum(y0(5:19));
% y0(5) = y0(5);

%% Start the simulation
y = y0; T = t_end;
t_cri_lst = unique(sort([tk_pre_lst, tk_post_lst]));
[t_lst, y_lst] = Simu_ode_scheme_multistep(0, T, t_cri_lst, 1, 1e-4, y, paras, cal_source);

c_MP_lst = y_lst(:, 6:19) * paras.MP_lst';
c_M_lst = y_lst(:, 6:19) * paras.M_lst';
c_E_lst= y_lst(:, 1);

% y_ss = y_lst(end, :);
% CaMKII_SS = y_ss(6:19) * paras.MP_lst' / (paras.n * paras.c0) * 100;
% E_SS = y_ss(1); 
% AMPA_SS = y_ss(20) / y0(20);

fprintf('%d finished! \n', k)
end

function [c_Ca] = Ca_func(sol, t)
y = deval(sol, t);
c_Ca = y(16);
end
