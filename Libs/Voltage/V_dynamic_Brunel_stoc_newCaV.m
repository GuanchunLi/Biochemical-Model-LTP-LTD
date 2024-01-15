function [dy_cir, J_Ca, J_CaV12, J_CaV13, J_NMDA] = V_dynamic_Brunel_stoc_newCaV(y_cir, x_AMPA, x_NMDA, I_stim, NMDA_factor, CaV_factor)

ratio_Cm = 1; % Normalization factor for capacitance

CaV12_factor = CaV_factor;
CaV13_factor = CaV_factor;

% Unpack the variables
Vs = y_cir(1); % mV Membrane voltage of spine head
s_AMPA = y_cir(2); % 1 Value of AMPA gating variable s
s_NMDA = y_cir(3); % 1 Value of NMDA gating variable s

y_CaV12 = y_cir(4:9);  % Variables of CaV12 model
y_CaV13 = y_cir(10:15); % Variables of CaV13 model

c_Ca = y_cir(16); % muM Intracelluar Calcium concentration
c_Ca_SS = y_cir(17); % muM Subspace Calcium concentration

s_Na = y_cir(18:19); % 1 Value of Sodium Current m, h
s_K = y_cir(20); % 1 Value of Potassium Current n

% Parameter setting (MAYBE NEEDS FINE TUNING!)
% cs = 0.05;    % muM Intracelluar Calcium concentration
g_CaV12 = 40; %  40, 182 mmol/(cm C) Strength of Ca current flux due to CaV12
g_CaV13 = 5;  %  100, 20 mmol/(cm C) Strength of Ca current flux
g_CaV12_SS = 0; %  mmol/(cm C) Strength of Ca current flux due to CaV12
g_CaV13_SS = 0;  %  mmol/(cm C) Strength of Ca current fluxdue to CaV13
F = 96.5; % C/mmol Faraday constant
vi = 1;   % (mu m)^3 Volume of spine head
alpha = F * vi * 1e-6; % pC / muM Constant (transforming flux to current)

g_Na = 0.7 * 1e3 / ratio_Cm; % nS Maximum sodium conductance
g_K = 1.3 * 1e3 / ratio_Cm; % nS Maximum potassium conductance
E_Na = 60; % mV Sodium reversal potential
E_K = -80; % mV Potassium reversal potental

beta_NMDA = 1e-3; % 1 Buffer factor for NMDA Ca current
beta_CaV = 1e-1 / 10; % 1 Buffer factor for CaV Ca current

% Extra inputs for the rates
I_AMPA = AMPA_current(s_AMPA, Vs);
[I_NMDA, I_NMDA_Ca] = NMDA_current(s_NMDA, Vs);
I_NMDA = I_NMDA * NMDA_factor;
I_NMDA_Ca = I_NMDA_Ca * NMDA_factor;
I_Na = - g_Na * s_Na(1).^3 * s_Na(2) * (Vs - E_Na); % pA Sodium current
I_K = - g_K * s_K.^4 * (Vs - E_K); % pA Potassium current
I_stim = I_stim / ratio_Cm + I_Na + I_K;

% CALCIUM LEFT FOR FUTURE WORK
Po_CaV12 = y_CaV12(1).^2;
Po_CaV13 = y_CaV13(1).^2;
J_Ca_12 = CaV_flux(Po_CaV12, g_CaV12, c_Ca, Vs) * CaV12_factor;
J_Ca_13 = CaV_flux(Po_CaV13, g_CaV13, c_Ca, Vs) * CaV13_factor;
J_Ca_12_SS = CaV_flux(Po_CaV12, g_CaV12_SS, c_Ca, Vs) * CaV12_factor;
J_Ca_13_SS = CaV_flux(Po_CaV13, g_CaV13_SS, c_Ca, Vs) * CaV13_factor;

J_Ca_NMDA = I_NMDA_Ca / (2*alpha);

J_Ca = (J_Ca_12 + J_Ca_13) * beta_CaV + J_Ca_NMDA * beta_NMDA;
J_CaV12 = J_Ca_12 * beta_CaV;
J_CaV13 = J_Ca_13 * beta_CaV;
J_NMDA = J_Ca_NMDA * beta_NMDA;
J_Ca_SS = (J_Ca_12_SS + J_Ca_13_SS) * beta_CaV + J_Ca_NMDA * beta_NMDA;
J_ER = 0; % Flux from ER release

% Transfrom the Ca flux to Ca current
I_CaV = 2*alpha*(J_Ca_12 + J_Ca_13); 

% Compute the rates
dVs = V_dynamic(Vs, I_AMPA, I_NMDA, I_CaV, I_stim);

ds_AMPA = AMPA(s_AMPA, x_AMPA);
ds_NMDA = NMDA(s_NMDA, x_NMDA);

dy_CaV12 = CaV12_New(y_CaV12, Vs);
dy_CaV13 = CaV13_New(y_CaV13, Vs);

dCa = cs_dynamic(c_Ca, J_Ca);
dCa_SS = cp_dynamic(c_Ca_SS, c_Ca, J_Ca_SS, J_ER);

ds_Na = Na_channel(s_Na, Vs);
ds_K = K_channel(s_K, Vs);
% Pack the rates
dy_cir = [dVs; ds_AMPA; ds_NMDA; dy_CaV12; dy_CaV13; dCa; dCa_SS; ds_Na; ds_K];
end