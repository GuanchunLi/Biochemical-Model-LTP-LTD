%% Initialization & Parameter setting
clear;

% addpath([pwd, '/Libs/']);
addpath([pwd, '/../Parameters']);

% Set up all parameters for the model
paras.n = 6;    % Number of monomers in each CaMKII (consider only one ring)
paras.v = 1;    % Reaction speed scale (reaction_rate = base_rate * v)
paras = Param_CaM(paras);
paras = Param_Concent(paras);
paras = Param_FastEqual(paras);
paras = Param_ReactRate(paras);
paras = Param_PhosTrans(paras);
paras = Param_PP1(paras);
paras.display = 0;

% Define the time course of calcium source
% cal_source = @(t) 0.0994;
freq_cal = 10; t_induce = 100; duration = 2;
cal_source = @(t) 0.055 + 0.0 * sum(exp(-(t-t_induce-(1:(min(t-t_induce, duration)*freq_cal))/freq_cal)/0.1));

paras = Param_RateModel(paras);

p_E = 0.1;
c_CaM4_tot = 0.2189;
% c_CaM4_tot = 0.0981;

% y_CaMKII = [11.9863302632848;8.71321895119819e-05;0.000630393027086376;0.000889218818508665;0.000763800605366619;0.000380250342967610;0.00125056350499765;0.00107746932660862;0.00107765584899687;0.000309953717241625;0.00175905645891763;0.00152223496408135;0.000763299625241066;0.00248471525829757;0.000673991839194866];

PP1_func = @(t, y) PP1_dynamic(t, y, p_E, c_CaM4_tot, paras);

y0 = [20; 5.0239];
T = 10;
[t_lst, y_lst] = ode45(PP1_func, [0, T], y0);
figure(1);
plot(t_lst, y_lst(:, 1));

c_E0_lst = 0:1e-2:50;
d_c_E0_lst = 0*c_E0_lst;
for i = 1:length(c_E0_lst)
    c_E0_i = c_E0_lst(i);
    dy_i = PP1_func(10, [c_E0_i; 5]);
    d_c_E0_i = dy_i(1);
    d_c_E0_lst(i) = d_c_E0_i;
end
figure(2);
plot(c_E0_lst, d_c_E0_lst);

%%
figure(3);
K_E0 = 1;
p_E = 1;
c_E0_lst = 0:1e-2:1;
paras.c_I_base = 0.1;
k_c_E = 10000;
K_c_E = 30;
c_I = y_lst(end, 2);
plot(c_E0_lst, paras.alpha_IE * (c_I - paras.c_I_base) * p_E * (c_E0_lst - paras.c_E0_base)./(K_E0 + c_E0_lst - paras.c_E0_base) );
hold on;
plot(c_E0_lst, (paras.beta_IE + k_c_E*(c_E0_lst- paras.c_E0_base)./(K_c_E + paras.c_E_tot - c_E0_lst)).*(paras.c_E_tot - c_E0_lst));