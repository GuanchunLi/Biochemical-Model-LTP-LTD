function [R_CaV] = Rate_CaV13_channel(n_CaV, V, cp)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% n_CaV -- number of CaV13 channels in different states: 7*1 array
% R_CaV-- reaction rates of CaV13 channels: 20*1 array

% Parameter setting
tau_po = 1;  % ms  Time constant of activation
kp0 = 3;    % muM  Threshold for Ca-induced inactivation
cp_b = 3;   % muM  Threshold for Ca dependence of transition rate k6

r1 = 0.3;  % 1/ms  Opening rate
r2 = 3;    % 1/ms  Closing rate
s1_p = 0.00195;  % 1/ms Inactivation rate
k1_p = 0.00413;  % 1/ms Inactivation rate
k2 = 0.0001;     % 1/ms Inactivation rate
k2_p = 0.00224;  % 1/ms Inactivation rate
T_Ba = 450;      % ms  Time constant

% Reaction rates
V_dis = 30;  % mV displacement of voltage from CaV1.2 model

po_inf = 1 ./ (1+exp(-(V+V_dis)/8));
alpha = po_inf / tau_po;
beta = (1 - po_inf) / tau_po;

f_cp = 1 / (1 + (kp0/cp)^3);
s1 = 0.02 * f_cp;
k1 = 0.03 * f_cp;
s2 = s1 * (k2/k1) * (r1/r2);
s2_p = s1_p * (k2_p/k1_p) * (r1/r2);

k3 = exp(-(V+V_dis+40)/3) ./ (3 * (1+exp(-(V+V_dis+40)/3)));
k3_p = k3;

R_V = 10 + 4954 * exp((V+V_dis)/15.6);
Pr = 1 - 1 / (1+exp(-(V+V_dis+40)/4));
Ps = 1 / (1+exp(-(V+V_dis+40)/11.32));
T_Ca = (78.0329 + 0.1*(1+cp/cp_b)^4) / (1 + (cp/cp_b)^4);
tau_Ca = (R_V - T_Ca) * Pr + T_Ca;
tau_Ba = (R_V - T_Ba) * Pr + T_Ba;

k5 = (1 - Ps) / tau_Ca;
k6 = f_cp * Ps / tau_Ca;
k5_p = (1 - Ps) / tau_Ba;
k6_p = Ps / tau_Ba;

k4 = k3 * (alpha/beta) * (k1/k2) * (k5/k6);
k4_p = k3_p * (alpha/beta) * (k1_p/k2_p) * (k5_p/k6_p);

R_lst = [s1, s2, r1, r2, s1_p, s2_p, ...
         k1, k2, k1_p, k2_p, ...
         k3, k4, alpha, beta, k3_p, k4_p, ...
         k5, k6, k5_p, k6_p]' * 1e3; % 1/s
n_lst = [n_CaV(1), n_CaV(2), n_CaV(3), n_CaV(1), n_CaV(1), n_CaV(4), ...
         n_CaV(3), n_CaV(2), n_CaV(3), n_CaV(4), ...
         n_CaV(2), n_CaV(5), n_CaV(6), n_CaV(3), n_CaV(4), n_CaV(7), ...
         n_CaV(5), n_CaV(6), n_CaV(7), n_CaV(6)]';

R_CaV = n_lst .* R_lst;

end

