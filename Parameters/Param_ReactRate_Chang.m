function [paras] = Param_ReactRate_Chang(paras)
%PARAM_REACTRATE Set parameters for the reaction rates
v = paras.v;

% Reaction rate for the self phosphorylation (neigh-phos)
paras.alpha_A = 0*v; paras.beta_A = 0*v;
paras.alpha_AC = 0*v; paras.beta_AC = 0*v;
paras.alpha_E = 200*v; paras.beta_E = 0*v;
% Reaction rate for the self phosphorylation (non-neigh-phos)
% paras.alpha_A = 0.001*v; paras.beta_A = 0*v;
% paras.alpha_AC = 2*v; paras.beta_AC = 0*v;
% paras.alpha_E = 3*v; paras.beta_E = 0*v;

% Reaction rate for the neighbouring phosphorylation
paras.alpha_AO = 100*v; paras.beta_AO = 0*v; % Cat-subunit is unphosphorylated
paras.alpha_AP = 100*v; paras.beta_AP = 0*v; % Cat-subunit is phosphorylated

% Reaction rate for the open and close (neigh-phos)
v2 = 10*v;
paras.mu = 1.1e-2*v2; paras.nu = (5e2)*v2;
% Reaction rate for the open and close (non-neigh-phos)
% paras.mu = 1.1e-2*v; paras.nu = (5*10^12)*v;
end

