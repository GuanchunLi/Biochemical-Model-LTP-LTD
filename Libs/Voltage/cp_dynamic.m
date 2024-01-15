function [d_cp] = cp_dynamic(cp, cs, J_Ca, J_ER)
% Dynamics of calcium concentration cp (timescale is 1s)
% Input:     cp -- muM calcium concentration
%          J_Ca -- muM/s Flux of Calcium through CaV1 channel
%          J_ER -- muM/s Flux through ER release
% Output:  d_cp -- muM/s Change rate of calcium concentration cp

% Parameter setting
% cs = 0.05;  % muM
tau_s = 0.5; % 0.5;  % ms Dyadic junction-submembrane diffusion time constant

d_cp = J_Ca + J_ER - (cp - cs) * 10^3 / tau_s ; % muM/s
end