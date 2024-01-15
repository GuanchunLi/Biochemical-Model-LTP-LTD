function [dy_PP1] = PP1_dynamic_mod5(t, y_PP1, c_Ca, p_E, p_EP, dE_CaM, c_CaM4_tot, paras)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Unpack the variables
c_E0 = y_PP1(1); c_E = p_E * c_E0;
c_EP0 = y_PP1(3); c_EP = p_EP * c_EP0;
c_I = y_PP1(4); c_IE = y_PP1(2);

[c_PP1, c_PP1_Cdk5, c_PP1_PP1P] = PP_fast_solver(c_E0, c_EP0, c_Ca, paras);

dy_PP1 = 0 * y_PP1;
dy_PP1(1) = - dy_PP1(2) + paras.alpha_P_auto * c_PP1_PP1P ...
          - paras.alpha_Cdk5 * c_PP1_Cdk5 + dE_CaM + paras.d_E_base;
dy_PP1(3) = - (dy_PP1(1) + dy_PP1(2));
end

