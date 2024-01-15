function [dy] = CaV13_New(y, V)

% Parameter setting
m = y(1);
% Reaction rates
a = 1.6 / (1 + exp(-0.072*(V+16)));
b = - 0.02 * (V+39.4) / (1 - exp((V+39.4)/5.36));

m_inf = 1 / (1 + exp(-(V+39.4)/7));
tau_m = 1 / (a + b);

% Compute the rates
dy = 0*y;

dy(1) = (m_inf - m) / tau_m;
% dy(2) = alpha*C2 + k2*I1C + k2_p*I1B + r2*Po - (r1+beta+k1+k1_p)*C1;
% dy(3) = k1*C1 + k4*I2C + s1*Po - (k2+k3+s2)*I1C;
% dy(4) = k3*I1C + k6*C2 - (k4+k5)*I2C;
% dy(5) = k1_p*C1 + k4_p*I2B + s1_p*Po - (k2_p+k3_p+s2_p)*I1B;
% dy(6) = k3_p*I1B + k6_p*C2 - (k5_p+k4_p)*I2B;

dy = dy * 10^3; % Convert unit of dy to 1/s
end