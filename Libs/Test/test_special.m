%% Initialization & Parameter setting
clear;

addpath([pwd, '/../Parameters']);

% Set up all parameters for the model
paras.n = 6;    % Number of monomers in each CaMKII (consider only one ring)
paras.v = 1.3;    % Reaction speed scale (reaction_rate = base_rate * v)
paras = Param_CaM(paras);
paras = Param_Concent(paras);
paras = Param_FastEqual(paras);
paras = Param_ReactRate(paras);
paras = Param_PhosTrans(paras);
paras.display = 0;

% Define the time course of calcium source
cal_source = @(t) 0.1;

% Initialization for CaMKII (ODE version)
y0 = zeros(8, 1);
y0(8, :) = 3*0.9860;
y0(7, :) = 3*0.0519;
paras.MP_lst = [0, 1, 2, 3, 4, 5, 6];
paras.M_lst = paras.n - paras.MP_lst;

c0 = 3*1.0379;

paras.c0 = c0;

c_Ca = 0.1;
%% Test
c_Ms = 8; c_MPs = 5; c_Ca = 0.5;
% K_Ca = equal_constant_Ca(c_Ca, paras.K_1, paras.K_2, paras.K_3, paras.K_4);
% X1 = solve_fast_equal_special(paras.c_CaM, paras.c_E0, c_Ms, c_MPs, ...
%     paras.c_ATP, paras.c_ADP, paras.c_P, K_Ca, paras);
% X2 = solve_fast_equal(paras.c_CaM, paras.c_E0, c_Ms, c_MPs,...
%                      paras.c_ATP, paras.c_ADP, paras.c_P, K_Ca, paras);
% [p_M, p_M_ATP, p_M_CaM, p_MA_CaM, p_M_EP, ...
%           p_MP, p_MP_ADP, p_MP_CaM, p_MPA_CaM, p_MP_E] ...
% = compute_Prob(c_Ms, c_MPs, c_Ca, paras);

[X, flag] = solve_SS_non_phos(paras, c_Ca);