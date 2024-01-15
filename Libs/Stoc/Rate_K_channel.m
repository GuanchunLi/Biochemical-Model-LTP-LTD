function [R_K] = Rate_K_channel(n_K, V)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% n_K -- number of potassium channels in different states: 5*1 array
% R_K -- reaction rates of sodium channels: 8*1 array

K_inf = 1 ./ (1 + exp(-(V+30)/25)); % 1 n_inf
tau_K = 0.01 + 2.5 ./ (exp((V+30)/40) + exp(-(V+30)/50)); % ms

alpha_K = K_inf / tau_K * 1e3; % 1/s open rate
beta_K = (1 - K_inf) / tau_K * 1e3; % 1/s closing rate

R_lst = [4*alpha_K, beta_K, 3*alpha_K, 2*beta_K, 2*alpha_K, 3*beta_K, alpha_K, 4*beta_K]'; % 1/s
n_lst = [n_K(1), n_K(2), n_K(2), n_K(3), n_K(3), n_K(4), n_K(4), n_K(5)]';
R_K = n_lst .* R_lst;

end

