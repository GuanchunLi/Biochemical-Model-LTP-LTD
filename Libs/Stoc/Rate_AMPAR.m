function [R_AMPAR] = Rate_AMPAR(n_AMPAR)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% n_AMPAR -- number of AMPARs in different states: 4*1 array
% R_AMPAR -- reaction rates of AMPARs: 5*1 array

% Set the parameters
tau_a = 0.05;  % ms Time constant for AMPA gating variable x_AMPA
tau_s = 2;  % ms Time constant for AMPA gating variable s_AMPA
alpha_s = 1;   % 1/ms Effect magnitude of presynaptic variable x_AMPA

R_lst = [1/tau_a, 1/tau_a, 1/tau_s, 1/tau_s, alpha_s]' * 1e3; % 1/s
n_lst = [n_AMPAR(2), n_AMPAR(4), n_AMPAR(3), n_AMPAR(4), n_AMPAR(2)]';
R_AMPAR = n_lst .* R_lst;

end

