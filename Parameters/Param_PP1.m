function [paras] = Param_PP1(paras)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% PKA parameters
paras.k_PKA_0 = 0.00359;  % 1/s PKA base activity
paras.k_PKA = 100;  % 1/s Maximum CaM dependent PKA activity
paras.K_PKA = 0.11; % muM PKA half activity concentration
paras.n_PKA = 8; % 1 PKA Hill coefficient

% CaN parameters
paras.k_CaN_0 = 0.1;  % 1/s CaN base activity
paras.k_CaN = 20; % 1/s Maximum CaM dependent CaN activity
paras.K_CaN = 0.053; % muM CaN half activity concentration
paras.n_CaN = 3; % 1 CaN Hill coefficient

% I1P, PP1 reaction rate
paras.alpha_IE = 0; % 1/(s*muM) I1P, PP1 association rate
paras.beta_IE = 0; % 1/s I1P, PP1 dissociation rate

% PP1 auto-dephosphorylation parameters
paras.alpha_P_auto = 1e-2; % 1000 1/s PP1 auto-dephosphorylation rate
paras.alpha_P = 0; % 1/s PP1 dephosphorylation rate through Ca pathway

paras.alpha_Cdk5 = 10; % 1/s Cdk5-PP1 phosphorylation rate
paras.K_m50_EP = 100; % muM PP1 auto-dephosphorylation Michalis Menton Coef
paras.K_m50_E = 20; % muM Cdk5-PP1 phosphorylation Michalis Menton Coef

paras.d_E_base = 1e-2;  % muM/s Baseline dephosphorylation rate should be 1e-4

% Concentration rate
paras.c_I0 = 1; % muM Total I1 concentration

% Baseline activity of c_E0
paras.r_P_base = 0; % 1/s 

% paras.c_E0_base = 0.9; % muM

% Total concentration of Cdk5
paras.Cdk5_tot = 2; % muM

end

