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

color_lst = [[0.215686274509804 0.494117647058824 0.72156862745098];
    [0.894117647058824 0.101960784313725 0.109803921568627];
    [0.909803921568627 0.623529411764706 0.290196078431373];
    [0.56078431372549 0 0.980392156862745]];
%% Plot Calcium Dynamics
t_end = 40;
tk_pre_center = [15.1];
tk_pre_lst = [15:0.1:15.2]; % Pre-spike time list, unit: s
% tk_pre_center = [15.5];
% tk_pre_lst = [15:0.1:16]; % Pre-spike time list, unit: s

% dt_lst = [3000];
dc_lst = [0, 100];
figure('OuterPosition',[623 13 594.666666666667 918]);
axes1 = axes;
hold(axes1,'on');
hold on;
for k = 1:length(dc_lst)
    dc = dc_lst(k);
    [t_Ca_lst, Ca_lst, Ca_sys_lst] = run_simulation_k(tk_pre_lst, tk_pre_center, t_end, dc_lst(k), k, paras);
%     save(sprintf('Result/Result_04_28/BTSP/BTSP_dynamic_dt=%dms.mat', dt), 'dt', ...
%      't_Ca_lst', 'Ca_lst', 't_lst', 'y_lst', 'c_MP_lst', 'c_E_lst');
    % fig = figure(k);
    plot(t_Ca_lst, Ca_lst, 'DisplayName', sprintf('%d pA', dc), 'LineWidth',2.5,...
    'Color', color_lst(k,:));
%     xlabel('Time (s) ', 'FontSize', 18);
%     ylabel('Calcium Concentration (\muM)', 'FontSize', 18);
%     saveas(fig, sprintf('Result/Figures/BTSP/BTSP_dynamic_dt=%dms.png', dt));
end

% 创建 ylabel
ylabel('Ca (\muM)');

% 创建 xlabel
xlabel('Time (s) ');

% 取消以下行的注释以保留坐标区的 X 范围
xlim(axes1,[12 20]);
% 取消以下行的注释以保留坐标区的 Y 范围
% ylim(axes1,[0.0014359752197185 1.40143597521972]);
hold(axes1,'off');
% 设置其余坐标区属性
set(axes1,'FontSize',28,'LineWidth',2.5,'XTick',...
    [14.45 15.45 16.45 17.45 18.45 19.45 20.45 21.45],'XTickLabel',...
    {'-1','0','1','2','3','4','5','6'},'XTickLabelRotation',45,'YTick',...
    [0 0.25 0.5 1 1.25],'YTickLabel',{'0','0.25','0.5','1','1.25'});
% 创建 legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.512068965517241 0.555533266351566 0.329310344827586 0.208700882117081],...
    'EdgeColor','none');

function [t_Ca_lst, Ca_lst, Ca_sys_lst] = run_simulation_k(tk_pre_lst, tk_pre_center, t_end, dc, k, paras)
tk_post_lst = max(tk_pre_lst) + 500/1000;
% tk_post_lst = [];

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
I_stim_func = @(t) 4000 * any(((tk_post_lst + 300/1e3) > t) .* (tk_post_lst < (t))) + dc;

circuit_model = @(t, y) V_dynamic_Brunel_stoc(y, x_AMPA_func(t),...
                           x_NMDA_func(t), I_stim_func(t), ... % 0.418, 2);
                           0.412, 2.02);

opts = odeset('RelTol',1e-7, 'AbsTol',1e-6, 'MaxStep',1e-3);
dc_model = @(t, y) V_dynamic_Brunel_stoc(y, 0, 0, I_stim_func(t), 0.41, 2);
[t_dc, y_dc] = ode15s(dc_model, [0, 5], paras.y0_Ca, opts);
y0 = paras.y0_Ca;
y0 = y_dc(end,:)';

sol_Ca = ode15s(circuit_model, [0, t_end], y0, opts);
t = sol_Ca.x; y = sol_Ca.y';
cal_source = @(t) Ca_func(sol_Ca, t);
t_Ca_lst = t; Ca_lst = y(:, 16); V_lst = y(:, 1); Ca_sys_lst = y;
end