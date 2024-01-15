function [R_NMDAR] = Rate_NMDAR(n_NMDAR)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% n_NMDAR -- number of NMDARs in different states: 4*1 array
% R_NMDAR -- reaction rates of NMDARs: 5*1 array

% Set the parameters
tau_a = 2;  % ms Time constant for NMDA gating variable x_NMDA
tau_s = 80;  % ms Time constant for NMDA gating variable s_NMDA
alpha_s = 1;   % 1/ms Effect magnitude of presynaptic variable x_NMDA

R_lst = [1/tau_a, 1/tau_a, 1/tau_s, 1/tau_s, alpha_s]' * 1e3; % 1/s
n_lst = [n_NMDAR(2), n_NMDAR(4), n_NMDAR(3), n_NMDAR(4), n_NMDAR(2)]';
R_NMDAR = n_lst .* R_lst;

end

