function [dy_cir, J_Ca] = circuit_dynamic(y_cir, x_AMPA, x_NMDA, cp, I_stim)

% Unpack the variables
Vs = y_cir(1); % mV Membrane voltage of spine head
Vd = y_cir(2); % mV Membrane voltage of dendrite
s_AMPA = y_cir(3); % 1 Value of AMPA gating variable s
s_NMDA = y_cir(4); % 1 Value of NMDA gating variable s

y_CaV12 = y_cir(5:10);  % Variables of CaV12 model
y_CaV13 = y_cir(11:16); % Variables of CaV13 model

% Parameter setting (MAYBE NEEDS FINE TUNING!)
R_N = 500; % M_Omega Resistance of spine neck
g_LD = 2.6e-3; % muS Leak conductance of dendrite
Ed = -70; % mV  Baseline potential of dendrite (from soma?)
C_d = 0.052;% nF Dendrite capacitance
cs = 0;    % muM Intracelluar Calcium concentraion
g_CaV12 = 182*50; %  mmol/(cm C) Strength of Ca current flux due to CaV12
g_CaV13 = 20*50;  %  mmol/(cm C) Strength of Ca current flux due to CaV13
F = 96.5; % C/mmol Faraday constant
vi = 0.8;   % (mu m)^3 Volume of spine head
alpha = F * vi; % Constant (transforming flux to current)

% Extra inputs for the rates
I_AMPA = AMPA_current(s_AMPA, Vs);
[I_NMDA, I_NMDA_Ca] = NMDA_current(s_NMDA, Vs);
I_N = (Vs - Vd) / R_N * 10^3;  % pA Current through spine neck
I_d = (Vd - Ed) * g_LD * 10^3;  % pA Leakage current through dendrite
I_stim = I_stim - I_N;

% CALCIUM LEFT FOR FUTURE WORK
Po_CaV12 = 1 - sum(y_CaV12);
Po_CaV13 = 1 - sum(y_CaV13);
J_Ca_12 = CaV_flux(Po_CaV12, g_CaV12, cs, Vs);
J_Ca_13 = CaV_flux(Po_CaV13, g_CaV13, cs, Vs);

J_Ca_NMDA = I_NMDA * 10^6 / (2*alpha) * 0.1;
J_Ca = J_Ca_12 + J_Ca_13 + J_Ca_NMDA * 0;

% Transfrom the Ca flux to Ca current
I_CaV = 2*alpha*(J_Ca_12 + J_Ca_13)/(10^6); 

%I_CaV = 0;
%I_NMDA = 0;
%I_AMPA = 0;

% Compute the rates
dVs = V_dynamic(Vs, I_AMPA, I_NMDA, I_CaV, I_stim);
dVd = (I_N - I_d) / C_d;
ds_AMPA = AMPA(s_AMPA, x_AMPA);
ds_NMDA = NMDA(s_NMDA, x_NMDA);

dy_CaV12 = CaV12(y_CaV12, cp, Vs);
dy_CaV13 = CaV12(y_CaV13, cp, Vs);

% Pack the rates
dy_cir = [dVs; dVd; ds_AMPA; ds_NMDA; dy_CaV12; dy_CaV13];
end