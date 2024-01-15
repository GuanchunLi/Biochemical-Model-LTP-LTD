function [paras] = Param_PP1_mod4(paras)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% PP1 auto-dephosphorylation parameters
paras.alpha_P_auto = 6; % 1000 1/s PP1 auto-dephosphorylation rate
paras.K_PP1_PP1P = 1e-1; % muM  PP1-PP1P binding parameter

% Cdk5-PP1 phosphorylation parameters
paras.alpha_Cdk5 = 2; % 1/s Cdk5-PP1 phosphorylation rate
paras.K_PP1_Cdk5 = 1e-2; % muM PP1-Cdk5 binding parameter
paras.v0 = 1;

% Baseline activity of dephosphorylation
paras.d_E_base = 1;  % muM/s Baseline dephosphorylation rate

% Baseline activity of c_E0
paras.r_P_base = 0; % 1/s 

% paras.c_E0_base = 0.9; % muM

% Total concentration of Cdk5
paras.Cdk5_tot = 50; % muM

end

