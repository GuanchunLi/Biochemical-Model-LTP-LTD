function [dy] = CaV_model_Mahajan(y, Vs)
cp = y(1);
cs = y(2);
y_CaV12 = y(3:8);
Po_CaV12 = 1 - sum(y_CaV12);

CaV_factor = 2.02;
beta_CaV = 1e-1 / 10; % 1 Buffer factor for CaV Ca current

g_CaV12 = 82 * 10; %  40, 182 mmol/(cm C) Strength of Ca current flux due to CaV12
% g_CaV13 = 20 * 10;  %  100, 20 mmol/(cm C) Strength of Ca current flux
g_CaV12_SS = 50*g_CaV12; %  mmol/(cm C) Strength of Ca current flux due to CaV12
% g_CaV13_SS = 50*g_CaV13;  %  mmol/(cm C) Strength of Ca current fluxdue to CaV13

J_Ca_12 = CaV_flux(Po_CaV12, g_CaV12, cs, Vs) * CaV_factor;
% J_Ca_13 = CaV_flux(Po_CaV13, g_CaV13, c_Ca, Vs) * CaV_factor;
J_Ca_12_SS = CaV_flux(Po_CaV12, g_CaV12_SS, cs, Vs) * CaV_factor;
% J_Ca_13_SS = CaV_flux(Po_CaV13, g_CaV13_SS, c_Ca, Vs) * CaV_factor;

J_ER = 0;

J_Ca = J_Ca_12 * beta_CaV;
J_Ca_SS = J_Ca_12_SS * beta_CaV;

d_cs = cs_dynamic(cs, J_Ca);
d_cp = cp_dynamic(cp, cs, J_Ca_SS, J_ER);
dy_CaV12 = CaV12(y_CaV12, cp, Vs);

dy = [d_cp; d_cs; dy_CaV12];
end