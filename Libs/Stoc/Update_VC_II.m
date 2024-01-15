function [y_VC_new, dy_VC] = Update_VC_II(y_VC, n_channels, I_stim, dt, NMDA_factor, CaV_factor)

% Unpack the variables
Vs = y_VC(1);
c_Ca = y_VC(2);
c_Ca_SS = y_VC(3);
s_Na = y_VC(4:5);
s_K = y_VC(6);

CaV12_factor = CaV_factor;
CaV13_factor = CaV_factor;

% Unpack the variables
n_AMPAR = n_channels(1:4); % 1 
n_NMDAR = n_channels(5:8); % 1
n_CaV12 = n_channels(9:15); % 1
n_CaV13 = n_channels(16:22); % 1
% n_Na = n_channels(23:30); % 1
% n_K = n_channels(31:35); % 1

% Parameter setting (MAYBE NEEDS FINE TUNING!)
% cs = 0.05;    % muM Intracelluar Calcium concentration
g_CaV12 = 82*10; %  40, 182 mmol/(cm C) Strength of Ca current flux due to CaV12
g_CaV13 = 20*10;  %  100, 20 mmol/(cm C) Strength of Ca current flux
g_CaV12_SS = 50*g_CaV12; %  mmol/(cm C) Strength of Ca current flux due to CaV12
g_CaV13_SS = 50*g_CaV13;  %  mmol/(cm C) Strength of Ca current fluxdue to CaV13
F = 96.5; % C/mmol Faraday constant
vi = 1;   % (mu m)^3 Volume of spine head
alpha = F * vi * 1e-6; % pC / muM Constant (transforming flux to current)

g_Na = 0.7 * 1e3; % nS Maximum sodium conductance
g_K = 1.3 * 1e3; % nS Maximum potassium conductance
E_Na = 60; % mV Sodium reversal potential
E_K = -80; % mV Potassium reversal potental

beta_NMDA = 1e-3; % 1 Buffer factor for NMDA Ca current
beta_CaV = 1e-1 / 10; % 1 Buffer factor for CaV Ca current

% Translate the numbers to the open probability
s_AMPA = (n_AMPAR(3) + n_AMPAR(4)) / sum(n_AMPAR);
s_NMDA = (n_NMDAR(3) + n_NMDAR(4)) / sum(n_NMDAR);

% open probablity of Na & K channels
p_Na = s_Na(1).^3 * s_Na(2);
p_K = s_K.^4;
% p_Na = n_Na(8) / sum(n_Na);
% p_K = n_K(5) / sum(n_K);

% Extra inputs for the rates
I_AMPA = AMPA_current(s_AMPA, Vs);
[I_NMDA, I_NMDA_Ca] = NMDA_current(s_NMDA, Vs);
I_NMDA = I_NMDA * NMDA_factor;
I_NMDA_Ca = I_NMDA_Ca * NMDA_factor;
I_Na = - g_Na * p_Na * (Vs - E_Na); % pA Sodium current
I_K = - g_K * p_K * (Vs - E_K); % pA Potassium current
I_stim = I_stim + I_Na + I_K;

% CALCIUM LEFT FOR FUTURE WORK
Po_CaV12 = n_CaV12(1) / sum(n_CaV12);
Po_CaV13 = n_CaV13(1) / sum(n_CaV13);
J_Ca_12 = CaV_flux(Po_CaV12, g_CaV12, c_Ca, Vs) * CaV12_factor;
J_Ca_13 = CaV_flux(Po_CaV13, g_CaV13, c_Ca, Vs) * CaV13_factor;
J_Ca_12_SS = CaV_flux(Po_CaV12, g_CaV12_SS, c_Ca, Vs) * CaV12_factor;
J_Ca_13_SS = CaV_flux(Po_CaV13, g_CaV13_SS, c_Ca, Vs) * CaV13_factor;

J_Ca_NMDA = I_NMDA_Ca / (2*alpha);

J_Ca = (J_Ca_12 + J_Ca_13) * beta_CaV + J_Ca_NMDA * beta_NMDA;
% J_CaV12 = J_Ca_12 * beta_CaV;
% J_CaV13 = J_Ca_13 * beta_CaV;
% J_NMDA = J_Ca_NMDA * beta_NMDA;
J_Ca_SS = (J_Ca_12_SS + J_Ca_13_SS) * beta_CaV + J_Ca_NMDA * beta_NMDA;
J_ER = 0; % Flux from ER release

% Transfrom the Ca flux to Ca current
I_CaV = 2*alpha*(J_Ca_12 + J_Ca_13); 

% Compute the rates
dVs = V_dynamic(Vs, I_AMPA, I_NMDA, I_CaV, I_stim);
dCa = cs_dynamic(c_Ca, J_Ca);
dCa_SS = cp_dynamic(c_Ca_SS, c_Ca, J_Ca_SS, J_ER);
ds_Na = Na_channel(s_Na, Vs);
ds_K = K_channel(s_K, Vs);

dy_VC = 0 * y_VC;
dy_VC(1) = dVs;
dy_VC(2) = dCa;
dy_VC(3) = dCa_SS;
dy_VC(4:5) = ds_Na;
dy_VC(6) = ds_K;

% Vs_new = Vs + dVs * dt;
% c_Ca_new = c_Ca + dCa * dt;
% c_Ca_SS_new = c_Ca_SS + dCa_SS * dt;

% Pack the variables
y_VC_new = y_VC + dy_VC * dt;

end