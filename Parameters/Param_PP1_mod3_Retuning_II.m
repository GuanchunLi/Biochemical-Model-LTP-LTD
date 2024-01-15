function [paras] = Param_PP1_mod3_Retuning_II(paras)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% PP1 auto-dephosphorylation parameters
paras.alpha_P_auto = 0.5 * 2.95/5; % 1000 1/s PP1 auto-dephosphorylation rate
paras.alpha_P = 4 / 5 * 1.12; % 1/s PP1 dephosphorylation rate through Ca pathway

paras.K_m50_EP = 1; % muM PP1 auto-dephosphorylation Michalis Menton Coef

paras.d_E_base = 1e-4;  % 1e-4 muM/s Baseline dephosphorylation rate should be 1e-4

% Baseline activity of c_E0
paras.r_P_base = 0; % 1/s 

end

