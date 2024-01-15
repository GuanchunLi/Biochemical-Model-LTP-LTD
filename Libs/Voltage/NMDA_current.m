function [I_NMDA, I_NMDA_Ca] = NMDA_current(s_NMDA, V)
% Compute the current through NMDA channel
% Input: s_NMDA -- 1 Value of NMDA gating variable s
%         V     -- mV  Membrane potential
% Output: I_NMDA -- pA Total Current through NMDA channel
%      I_NMDA_Ca -- pA Calcium Current through NMDA channel

% The unit of I_NMDA(_Ca) is pA

% Set the parameters
g_NMDA = 4.5e-1 * 2;   % nS Maximum NMDA-R current conductance
E_NMDA = 0;        % mV  Total NMDA-mediated current reversal potentia
E_Ca = 140;   % mV  Calcium reversal potential
c_Mg = 1;          % mM  Concentration of magnesium

B_V = 1 / (1 + exp(-0.062*V) * c_Mg / 3.57); % Gating of Magnesium

% Compute the current
I_NMDA = -g_NMDA * s_NMDA * B_V * (V - E_NMDA); % pA total current through NMDA channel
I_NMDA_Ca = -g_NMDA * s_NMDA * B_V * (V - E_Ca); % pA Calcium current through NMDA channel

end