function [R_Na] = Rate_Na_channel(n_Na, V)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% n_Na -- number of sodium channels in different states: 8*1 array
% R_Na -- reaction rates of sodium channels: 20*1 array

m_inf = 1 ./ (1 + exp(-(V+30)/8.5)); % 1 m_inf
tau_m = 0.1; % ms time constant of m
alpha_m = m_inf / tau_m * 1e3; % 1/s open rate
beta_m = (1 - m_inf) / tau_m * 1e3; % 1/s closing rate

h_inf = 1 ./ (1 + exp((V+44.1)/7)); % 1 h_inf
tau_h = 1 + 3.5 ./ (exp((V+35)/4) + exp(-(V+35)/25)); % ms time constant of h
alpha_h = h_inf / tau_h * 1e3; % 1/s open rate
beta_h = (1 - h_inf) / tau_h * 1e3; % 1/s closing rate

R_lst = [3*alpha_m, beta_m, 2*alpha_m, 2*beta_m, alpha_m, 3*beta_m, ...
         3*alpha_m, beta_m, 2*alpha_m, 2*beta_m, alpha_m, 3*beta_m, ...
         alpha_h, beta_h, alpha_h, beta_h, ...
         alpha_h, beta_h, alpha_h, beta_h]'; % 1/s
n_lst = [n_Na(1), n_Na(2), n_Na(2), n_Na(3), n_Na(3), n_Na(4), ...
         n_Na(5), n_Na(6), n_Na(6), n_Na(7), n_Na(7), n_Na(8), ...
         n_Na(1), n_Na(5), n_Na(2), n_Na(6), ...
         n_Na(3), n_Na(7), n_Na(4), n_Na(8)]';

R_Na = n_lst .* R_lst;

end

