function [dy_PP1] = PP1_dynamic_mod1(t, y_PP1, c_Ca, p_E, p_EP, dE_CaM, c_CaM4_tot, paras)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Unpack the variables
c_E0 = y_PP1(1); c_E = p_E * c_E0;
c_EP0 = y_PP1(3); c_EP = p_EP * c_EP0;
c_I = y_PP1(4); c_IE = y_PP1(2);

% Compute the calcium dependent rates of PP1 dephosphrylation
% alpha_P_auto_Ca = alpha_P_auto_func(c_Ca, paras);
K_m50_EP = paras.K_m50_EP;
K_m50_E = paras.K_m50_E;
% alpha_P_Ca = alpha_P_func(c_Ca, paras);
Cdk5_a_Ca = Cdk5_a_func(c_Ca, paras);
% paras.c_I_base = 0.1; paras.dE_base = 1;
% d_E_base = paras.d_E_base * (1 + 0 * (c_Ca > 1.4));

dy_PP1 = 0 * y_PP1;
% dy_PP1(2) = paras.alpha_IE * c_I * c_E - paras.beta_IE * c_IE;
% dy_PP1(4) = - dy_PP1(2) - v_CaN * c_I + v_PKA * paras.c_I0;
dy_PP1(1) = - dy_PP1(2) + paras.alpha_P_auto * c_EP0 / (K_m50_EP + c_EP0) * c_E ...
          - c_E/(c_E + K_m50_E) * Cdk5_a_Ca * paras.alpha_Cdk5 + dE_CaM - paras.r_P_base*c_E + paras.d_E_base;
dy_PP1(3) = - (dy_PP1(1) + dy_PP1(2));
end

