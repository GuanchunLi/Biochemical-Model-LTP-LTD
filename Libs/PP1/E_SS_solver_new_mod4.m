function [c_MP] = E_SS_solver_new_mod4(y_E_SS, paras, c_Ca)
c_E0 = y_E_SS(1);
c_EP0 = y_E_SS(3);

solver_eq = @(c_E) E_SS_solver_eq(c_E, c_E0, c_EP0, c_Ca, paras);
if solver_eq(c_E0) > 0
   c_E = fzero(solver_eq, [0, c_E0]);
else
   c_E = c_E0 * (1-1e-10);
end
c_MP_E = c_E0 - c_E;
c_MP = paras.K_MP_E * c_MP_E / c_E;
end


function [d] = E_SS_solver_eq(c_E, c_E0, c_EP0, c_Ca, paras)
c_MP_E = c_E0 - c_E;
dE_CaM = - paras.alpha_E * c_MP_E;
[c_PP1, c_PP1_Cdk5, c_PP1_PP1P] = PP_fast_solver(c_E, c_EP0, c_Ca, paras);
d = paras.alpha_P_auto * c_PP1_PP1P ...
  - paras.alpha_Cdk5 * c_PP1_Cdk5 + dE_CaM + paras.d_E_base;
end