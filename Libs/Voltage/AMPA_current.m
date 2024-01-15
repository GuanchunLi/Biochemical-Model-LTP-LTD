function [I_AMPA] = AMPA_current(s_AMPA, V)
% Compute the current through AMPA channel
% Input: s_AMPA -- 1 Value of AMPA gating variable s
%         V     -- mV  Membrane potential
% Output: I_AMPA -- pA Current through AMPA channel

% The unit of I_AMPA is pA

% Set the parameters
g_AMPA = 19.5; % / paras.ratio_Cm * paras.ratio_AMPAR; % 1.2;   % nS Maximum AMPA current conductance
E_AMPA = 0;        % mV  AMPA-mediated current reversal potential

% Compute the current
I_AMPA = -g_AMPA * s_AMPA * (V - E_AMPA); % pA current through AMPA channel
end