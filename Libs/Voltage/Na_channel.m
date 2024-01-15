function [ds_Na] = Na_channel(s_Na, V)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
m = s_Na(1); h = s_Na(2);

m_inf = 1 ./ (1 + exp(-(V+30)/8.5)); % 1 m_inf
tau_m = 0.1; % ms time constant of m

h_inf = 1 ./ (1 + exp((V+44.1)/7)); % 1 h_inf
tau_h = 1 + 3.5 ./ (exp((V+35)/4) + exp(-(V+35)/25)); % ms time constant of h

dm = (m_inf - m) / tau_m * 1e3; % 1/s
dh = (h_inf - h) / tau_h * 1e3; % 1/s

ds_Na = [dm; dh];

end

