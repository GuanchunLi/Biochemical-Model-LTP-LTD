function [c_MP] = E_SS_solver_new_mod5(y_E_SS, paras, c_Ca)
c_E0 = y_E_SS(1);
c_EP0 = y_E_SS(3);

[c_PP1, c_PP1_Cdk5, c_PP1_PP1P] = PP_fast_solver(c_E0, c_EP0, c_Ca, paras);
dE_CaM = -(paras.alpha_P_auto * c_PP1_PP1P ...
  - paras.alpha_Cdk5 * c_PP1_Cdk5 + paras.d_E_base);
c_MP_E = - dE_CaM / paras.alpha_E;
c_E = c_E0 - c_MP_E;
c_MP = paras.K_MP_E * c_MP_E / c_E;
end
