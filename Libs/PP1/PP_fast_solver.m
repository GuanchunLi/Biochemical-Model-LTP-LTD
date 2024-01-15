function [c_PP1, c_PP1_Cdk5, c_PP1_PP1P] = PP_fast_solver(c_E, c_EP0, c_Ca, paras)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
v_Ca = Cdk5_Ca_func(c_Ca, paras);
c_PP1_func = @(c_PP1) c_PP1 * (1 + ...
    + (1/paras.K_PP1_Cdk5) * paras.Cdk5_tot ./ (1+v_Ca+c_PP1/paras.K_PP1_Cdk5) ...
    + (1/paras.K_PP1_PP1P) * c_EP0 / (1+c_PP1/paras.K_PP1_PP1P)) - c_E;
c_PP1 = fzero(c_PP1_func, [0, c_E]);
c_PP1_Cdk5 = (c_PP1/paras.K_PP1_Cdk5) * paras.Cdk5_tot ./ (1+v_Ca+c_PP1/paras.K_PP1_Cdk5);
c_PP1_PP1P = (c_PP1/paras.K_PP1_PP1P) * c_EP0 / (1+c_PP1/paras.K_PP1_PP1P);
end

