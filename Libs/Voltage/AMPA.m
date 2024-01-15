function [ds_AMPA] = AMPA(s_AMPA, x_AMPA)
% Dynamics of AMPA channel given x_AMPA trace (timescale is 1s)
% Input: s_AMPA -- 1 Value of AMPA gating variable s
%        x_AMPA -- 1 Value of AMPA gating variable x
% Output: ds_AMPA  -- 1/s Change rate of AMPA gating variable s

% Set the parameters
tau_AMPA = 2;  % ms Time constant for AMPA gating variable s_AMPA
alpha_s = 1;   % 1/s Effect magnitude of presynaptic variable x_AMPA

% Compute the rates
ds_AMPA = - s_AMPA / tau_AMPA + alpha_s * x_AMPA * (1 - s_AMPA); % 1 / ms

ds_AMPA = ds_AMPA * 10^3; % Convert unit of ds_AMPA to 1/s
end