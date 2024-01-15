function [dV] = V_dynamic(V, I_AMPA, I_NMDA, I_CaV, I_stim)
% Dynamics of membrane potential V (timescale is 1s)
% Input:      V -- mV Membrane potential
%        I_AMPA -- nA Current through AMPA channel
%        I_NMDA -- nA Total Current through NMDA channel
%         I_CaV -- nA Current through CaV1 channel
%        I_stim -- nA External input
% Output: dV  -- mV/s Change rate of membrane potential

ratio_Cm = 1; % paras.ratio_Cm;

% Set the parameters
g_L = 5 / ratio_Cm; %0.56e-3;   % nS  Maximum leak conductance
E_L = -70;  % mV Leak reversal potential making resting to be -70 mV (maybe change it later...)
C_m = 100 / ratio_Cm; % 1.13e-2;       % pF Whole cell capacitance
I_L = -g_L * (V - E_L); % pA Leakage current

% I_AMPA = AMPA_current(s_AMPA, V);
% [I_NMDA, I_NMDA_Ca] = NMDA_current(s_NMDA, V);

dV = (I_L + I_AMPA + I_NMDA + I_CaV + I_stim) / C_m; % mV/ms

dV = dV * 10^3; % Convert unit of dV to mV/s
end