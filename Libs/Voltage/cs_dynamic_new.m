function [d_cs] = cs_dynamic_new(cs, J_Ca)
% Dynamics of calcium concentration cs (timescale is 1s)
% Input:     cs -- muM calcium concentration
%          J_Ca -- muM/s Flux of Calcium
% Output:  d_cs -- muM/s Change rate of calcium concentration cp

% Parameter setting
c_Ca_base = 0.094; % 0.094;  % 0.06; 0.1 muM
tau_s = 12; % 0.5;  % ms Calcium decay time constant

d_cs = J_Ca - (cs - c_Ca_base) * 10^3 / tau_s ; % muM/s
end