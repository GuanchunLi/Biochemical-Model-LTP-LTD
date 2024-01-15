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
% paras = Param_ReactRate_Detailed(paras);
paras = Param_PhosTrans(paras);
paras = Param_PP1_mod3_Retuning(paras);
paras = Param_RateModel(paras);
paras = Param_Readout(paras);
paras = Param_InitVal_Retuning(paras);
paras.display = 0;
paras.ratio_Cm = 1; % 10; % 10;
paras.ratio_AMPAR = 1; % 0.0;
paras.ratio_NMDAR = 1; % 1.048; % 1.039;
paras.ratio_CaV = 1; % 1.01;

%% Generate Calcium source
t_end = 150;
% tk_pre_center = [15.45];
% tk_pre_lst = [15:0.1:15.9];
tk_pre_center = [15.45, 31.45];
tk_pre_lst = [15:0.1:15.9, 31:0.1:31.9]; % Pre-spike time list, unit: s

% parpool(32);

flag_lst = [1, 2, 3, 4];
% flag_lst = [2];
dt_lst = linspace(-6000, 6000, 32);

for flag = flag_lst
    N_trail = length(dt_lst);
    CaMKII_lst = zeros(N_trail, 1);
    E_lst = zeros(N_trail, 1);
    AMPA_lst = zeros(N_trail, 1);
    parfor k = 1:length(dt_lst)
        [CaMKII_SS, E_SS, AMPA_SS] = run_simulation_k(tk_pre_lst, tk_pre_center, t_end, dt_lst(k), k, paras, flag);
        CaMKII_lst(k) = CaMKII_SS;
        E_lst(k) = E_SS;
        AMPA_lst(k) = AMPA_SS;
    end

    fig1 = figure();
    plot(dt_lst, CaMKII_lst);
    hold on;
    plot(dt_lst, E_lst);
    % saveas(fig1, sprintf('BTSP_curve_Concent_flag=%d_Det.png', flag));
    
    fig2 = figure('OuterPosition',...
    [696.333333333333 182.333333333333 1112.66666666667 636.666666666667]);
    axes1 = axes;
    hold(axes1,'on');
    plot(dt_lst, AMPA_lst, '*-', 'LineWidth',3, 'Color',[0, 0, 0]);
    
    % 创建 ylabel
    ylabel('EPSP rel.');
    
    % 创建 xlabel
    xlabel('Timing Difference (ms)');
    
    hold(axes1,'off');
    % 设置其余坐标区属性
    set(axes1,'FontSize',28,'LineWidth',2.5);
    
    % save(sprintf('BTSP_curve_flag=%d_Det.mat', flag), 'dt_lst', 'CaMKII_lst', 'E_lst', 'AMPA_lst');
end

%% Aux functions

function [CaMKII_SS, E_SS, AMPA_SS] = run_simulation_k(tk_pre_lst, tk_pre_center, t_end, dt, k, paras, flag)
tk_post_lst = tk_pre_center - dt/1000;

sigma_NMDA = 3.3/100;
NMDA_factor_lst = 0 * tk_pre_lst + 1;
% NMDA_factor_lst = max(NMDA_factor_lst + sigma_NMDA * randn(1, length(tk_pre_lst)),0);
NMDA_factor = @(t) interp1([0,tk_pre_lst, t_end], [1,NMDA_factor_lst,1], t, 'previous');

sigma_CaV = 10/100;
CaV_factor_lst = 0 * tk_post_lst + 1;
% CaV_factor_lst = max(CaV_factor_lst + sigma_CaV * randn(1, length(tk_post_lst)),0);
CaV_factor = @(t) interp1([0,tk_post_lst, t_end], [1,CaV_factor_lst,1], t, 'previous');

tau_x_AMPA = 0.05/(10^3); alpha_x_AMPA = 1; 
x_AMPA_func = @(t) x_spike(t, tk_pre_lst, tau_x_AMPA, alpha_x_AMPA);
tau_x_NMDA = 2/(10^3); alpha_x_NMDA = 1;
x_NMDA_func = @(t) x_spike(t, tk_pre_lst, tau_x_NMDA, alpha_x_NMDA);
if flag == 1 % Control case
    I_stim_func = @(t) 2000 * any(((tk_post_lst + 300/1e3) > t) .* (tk_post_lst < (t))) + 42;
elseif flag == 2 % Depression of CaV1.3 for all times
    I_stim_func = @(t) 4000 * any((tk_post_lst < t-2/1e3) .* (tk_post_lst > t-302/1e3));
elseif flag == 3 % Depression of CaV1.3 only during ‘Priming’ (before plateau potential)
    tk_min = min(min(tk_pre_lst), min(tk_post_lst));
    I_stim_func = @(t) 4000 * any(((tk_post_lst + 300/1e3) > t) .* (tk_post_lst < (t))) + 40 * (t > tk_min);
elseif flag == 4 % Depression of CaV1.3 only after 'Priming' (after plateau potential)
    tk_max = max(max(tk_pre_lst), max(tk_post_lst));
    I_stim_func = @(t) 4000 * any(((tk_post_lst + 300/1e3) > t) .* (tk_post_lst < (t))) + 40  * (t < tk_max);
end
circuit_model = @(t, y) V_dynamic_Brunel_stoc(y, x_AMPA_func(t),...
                           x_NMDA_func(t), I_stim_func(t), ...
                           0.412 / paras.ratio_NMDAR , 2.02 * paras.ratio_CaV); 
                           % 0.358, 2.35);
                           % 0.41, 2.02);

y0 = paras.y0_Ca;

opts = odeset('RelTol',1e-7, 'AbsTol',1e-6, 'MaxStep',1e-3);
sol_Ca = ode15s(circuit_model, [0, t_end], y0, opts);
t = sol_Ca.x; y = sol_Ca.y';
cal_source = @(t) Ca_func(sol_Ca, t);

%% Initial state of CaMKII system

% Baseline SS
y0 = paras.y0_CaMKII;
paras.c0 = sum(y0(5:19));

%% Start the simulation
y = y0; T = t_end;
t_cri_lst = unique(sort([tk_pre_lst, tk_post_lst]));
[t_lst, y_lst] = Simu_ode_scheme_multistep(0, T, t_cri_lst, 1, 1e-3, y, paras, cal_source);

y_ss = y_lst(end, :);
CaMKII_SS = y_ss(6:19) * paras.MP_lst' / (paras.n * paras.c0) * 100;
E_SS = y_ss(1); 
AMPA_SS = y_ss(20) / y0(20);

fprintf('%d finished! \n', k)
end

function [c_Ca] = Ca_func(sol, t)
y = deval(sol, t);
c_Ca = y(16);
end