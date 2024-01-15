function [ds_NMDA] = NMDA(s_NMDA, x_NMDA)
% Dynamics of NMDA channel given x_NMDA trace (timescale is 1s)
% Input: s_NMDA -- 1 Value of NMDA gating variable s
%        x_NMDA -- 1 Value of NMDA gating variable x
% Output: ds_NMDA  -- 1/s Change rate of NMDA gating variable s

% Set the parameters
tau_NMDA = 80;  % ms Time constant for NMDA gating variable s_NMDA
alpha_s = 1;   % 1/ms Effect magnitude of presynaptic variable x_AMPA

% Compute the rates
ds_NMDA = - s_NMDA / tau_NMDA + alpha_s * x_NMDA * (1 - s_NMDA); % 1/ms

ds_NMDA = ds_NMDA * 10^3; % Convert unit of ds_NMDA to 1/s
end