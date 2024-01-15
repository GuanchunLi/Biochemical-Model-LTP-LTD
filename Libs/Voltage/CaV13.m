function [dy] = CaV13(y, cp, V)

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

% Unpack the variables
C2 = y(1);
C1 = y(2);
I1C = y(3);
I2C = y(4);
I1B = y(5);
I2B = y(6);
Po = 1 - (C1+C2+I1C+I2C+I1B+I2B);

% Compute the rates
dy = 0*y;

dy(1) = beta*C1 + k5*I2C + k5_p*I2B - (k6+k6_p+alpha)*C2;
dy(2) = alpha*C2 + k2*I1C + k2_p*I1B + r2*Po - (r1+beta+k1+k1_p)*C1;
dy(3) = k1*C1 + k4*I2C + s1*Po - (k2+k3+s2)*I1C;
dy(4) = k3*I1C + k6*C2 - (k4+k5)*I2C;
dy(5) = k1_p*C1 + k4_p*I2B + s1_p*Po - (k2_p+k3_p+s2_p)*I1B;
dy(6) = k3_p*I1B + k6_p*C2 - (k5_p+k4_p)*I2B;

dy = dy * 10^3; % Convert unit of dy to 1/s
end