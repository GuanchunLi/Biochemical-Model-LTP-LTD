function [J_Ca] = CaV_flux(Po_CaV, g_Ca, cs, V)
% Compute the flux of Calcium through CaV(1.2/1/3) channel
% Input: Po_CaV -- 1 Open probability of CaV1 channel
%         g_CaV -- mmol/(cm C) Strength of local Ca flux due to CaV1 channel
%         V     -- mV  Membrane potential
%         cs    -- muM Concentration of intercelluar calcium
% Output: J_Ca  -- muM/s Flux of Calcium through CaV1 channel

% Parameter setting
P_Ca = 0.00054; % cm/s Constant
% g_Ca_b = 9000; %  mmol/(cm C) Strength of local Ca flux due to CaV1 channels
F = 96.5; % C/mmol Faraday constant
R = 8.315; % J/(mol*K) Universal gas constant
T = 308;   % K Temperature
c_Ca_o = 1800; % muM External calcium concentration

a = V*F/(R*T); % unit-less

if V == 0
    i_Ca = (4*P_Ca*F^2/(R*T)) * (cs*exp(2*a) - 0.341*c_Ca_o) * (R*T/(2*F));
else
    i_Ca = (4*P_Ca*V*F^2/(R*T)) * (cs*exp(2*a) - 0.341*c_Ca_o) / (exp(2*a)-1);
end

J_Ca = -g_Ca * Po_CaV * i_Ca;
end