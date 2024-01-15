function [paras] = Param_PP1_mod6(paras)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% PP1 auto-dephosphorylation parameters
paras.alpha_P_auto = 1e-1; % 1000 1/s PP1 auto-dephosphorylation rate
paras.alpha_P = 4; % 1/s PP1 dephosphorylation rate through Ca pathway

paras.K_m50_EP = 0.1; % muM PP1 auto-dephosphorylation Michalis Menton Coef

paras.d_E_base = 1e-4;  % muM/s Baseline dephosphorylation rate should be 1e-4

% Baseline activity of c_E0
paras.r_P_base = 9.9e-2; % 1/s 

end