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
paras.ratio_Cm = 10; % 10;
paras.ratio_AMPAR = 1; % 0.6377;
paras.ratio_NMDAR = 1.048; % 1.048; % 1.0381; % 1.0451; % 1.039;
paras.ratio_CaV = 1; % 1.004; % 1.01;

%% Generate Calcium source
t_end = 120;
tk_pre_lst = [5:1:65]; % Pre-spike time list, unit: s

parpool(32);

dt_lst = linspace(-150, 150, 60);
% dt = 200;
N_dt = length(dt_lst);
N_trail = 128;
CaMKII_mat = zeros(N_trail, N_dt);
E_mat = zeros(N_trail, N_dt);
AMPA_mat = zeros(N_trail, N_dt);

for j = 1:N_dt
    dt = dt_lst(j);
    parfor k = 1:N_trail
        [CaMKII_SS, E_SS, AMPA_SS] = run_simulation_k(tk_pre_lst, t_end, dt, k, paras);
        CaMKII_mat(k, j) = CaMKII_SS;
        E_mat(k, j) = E_SS;
        AMPA_mat(k, j) = AMPA_SS;
    end
    fprintf('%d finished! \n', j);
end

fig1 = figure(1);
errorbar(dt_lst, mean(CaMKII_mat), std(CaMKII_mat));
hold on;
errorbar(dt_lst, mean(E_mat), std(E_mat));
saveas(fig1, 'STDP_curve_Concent_NewParam.png')
% 
fig2 = figure(2);
errorbar(dt_lst, mean(AMPA_mat-1), std(AMPA_mat-1));
saveas(fig2, 'STDP_curve_AMPA_NewParam.png')

save('STDP_curve_NewParam.mat', 'dt_lst', 'CaMKII_mat', 'E_mat', 'AMPA_mat');

%% Aux functions

function [CaMKII_SS, E_SS, AMPA_SS] = run_simulation_k(tk_pre_lst, t_end, dt, k, paras)
sigma_NMDA = 1/100; %3.3
NMDA_factor_lst = 0 * tk_pre_lst + 1;
NMDA_factor_lst = max(NMDA_factor_lst + sigma_NMDA * randn(1, length(tk_pre_lst)),0);
NMDA_factor = @(t) interp1([0,tk_pre_lst, t_end], [1,NMDA_factor_lst,1], t, 'previous');

tk_post_lst = tk_pre_lst - dt/1000;
sigma_CaV = 2/100; %10
CaV_factor_lst = 0 * tk_post_lst + 1;
CaV_factor_lst = max(CaV_factor_lst + sigma_CaV * randn(1, length(tk_post_lst)),0);
CaV_factor = @(t) interp1([0,tk_post_lst, max(tk_post_lst)+1e-1, t_end], [1,CaV_factor_lst,1, 1], t, 'previous');

tau_x_AMPA = 0.05/(10^3); alpha_x_AMPA = 1; 
x_AMPA_func = @(t) x_spike(t, tk_pre_lst, tau_x_AMPA, alpha_x_AMPA);
tau_x_NMDA = 2/(10^3); alpha_x_NMDA = 1;
x_NMDA_func = @(t) x_spike(t, tk_pre_lst, tau_x_NMDA, alpha_x_NMDA);
I_stim_func = @(t) 4000 * any((tk_post_lst > t+2/1e3) .* (tk_post_lst < t+3/1e3));

circuit_model = @(t, y) V_dynamic_Brunel_stoc(y, x_AMPA_func(t),...
                           x_NMDA_func(t), I_stim_func(t),... % , 2.02 * paras.ratio_CaV, paras
                           0.412*NMDA_factor(t)/paras.ratio_NMDAR, 2.02*CaV_factor(t)*paras.ratio_CaV, paras);
% 0.43 for stochatic
y0 = paras.y0_Ca;

opts = odeset('RelTol',1e-7, 'AbsTol',1e-6, 'MaxStep',1e-3);
% [t_dc, y_dc] = ode15s(@(t, y) V_dynamic_Brunel_stoc(y, 0, 0, 0, 0.412, 2.02), [0,3], paras.y0_Ca, opts);
% y0 = y_dc(end,:)';

sol_Ca = ode15s(circuit_model, [0, t_end], y0, opts);
t = sol_Ca.x; y = sol_Ca.y';
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