function [paras] = Param_ReactRate(paras)
%PARAM_REACTRATE Set parameters for the reaction rates
v = paras.v;

v = 1;
% Reaction rate for the self phosphorylation (neigh-phos)
paras.alpha_A = 0*v; paras.beta_A = 0*v;
paras.alpha_AC = 0*v; paras.beta_AC = 0*v;
paras.alpha_E = 10*v; paras.beta_E = 0*v;
% Reaction rate for the self phosphorylation (non-neigh-phos)
% paras.alpha_A = 0.001*v; paras.beta_A = 0*v;
% paras.alpha_AC = 2*v; paras.beta_AC = 0*v;
% paras.alpha_E = 3*v; paras.beta_E = 0*v;

% Reaction rate for the neighbouring phosphorylation
% paras.alpha_AO = 100*v; paras.beta_AO = 0*v; % Cat-subunit is unphosphorylated
% paras.alpha_AP = 100*v; paras.beta_AP = 0*v; % Cat-subunit is phosphorylated
paras.alpha_AO = 200*v; paras.beta_AO = 0*v; % Cat-subunit is unphosphorylated
paras.alpha_AP = 200*v; paras.beta_AP = 0*v; % Cat-subunit is phosphorylated

% Reaction rate for the open and close (neigh-phos)
v2 = 5;
paras.mu = 5e-1*v2; paras.nu = (5e4)*v2; % 5e4
% paras.mu = 1.1*v2; paras.nu = (5e6)*v2;
% Reaction rate for the open and close (non-neigh-phos)
% paras.mu = 1.1e-2*v; paras.nu = (5*10^12)*v;
end

