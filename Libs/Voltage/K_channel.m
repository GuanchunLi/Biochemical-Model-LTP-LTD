function [ds_K] = K_channel(s_K, V)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
n = s_K;

n_inf = 1 ./ (1 + exp(-(V+30)/25)); % 1 m_inf
tau_n = 0.01 + 2.5 ./ (exp((V+30)/40) + exp(-(V+30)/50)); % ms

dn = (n_inf - n) / tau_n * 1e3; % 1/s

ds_K = [dn];

end

