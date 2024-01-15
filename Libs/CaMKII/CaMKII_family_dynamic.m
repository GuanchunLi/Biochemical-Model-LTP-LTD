function [dy_CaMKII_family] = CaMKII_family_dynamic(t, y_CaMKII_family, c_Ca, paras)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Unpack the variables
y_PP1 = y_CaMKII_family(1:4);
y_CaMKII = y_CaMKII_family(5:19);
r_AMPA = y_CaMKII_family(20);
c_E0 = y_PP1(1);
c_IE = y_PP1(2);
c_EP0 = y_PP1(3);
% c_EP0 = paras.c_E_tot - c_E0 - c_IE;
c_CaM = paras.c_CaM;

r_MPs = dot(y_CaMKII(2:end), paras.MP_lst) / (paras.n*paras.c0);
r_PPs = c_E0 / paras.c_E_tot;

% Compute the reaction rates
[dy_CaMKII, c_lst, pM_lst, pMP_lst, p_E, p_EP, dE_CaM, c_CaM4_tot] ...
  = Reaction_ode_wi_transpp_RE(t, y_CaMKII, c_Ca, c_E0, c_EP0, c_CaM, paras);
[dy_PP1] = PP1_dynamic_mod3(t, y_PP1, c_Ca, p_E, p_EP, dE_CaM, c_CaM4_tot, paras);
[d_r_AMPA] = AMPA_phos(r_AMPA, r_MPs, r_PPs, paras);

% Pack the varabiles
dy_CaMKII_family = 0 * y_CaMKII_family;
dy_CaMKII_family(1:4) = dy_PP1;
dy_CaMKII_family(5:19) = dy_CaMKII;
dy_CaMKII_family(20) = d_r_AMPA;

dy_CaMKII_family = dy_CaMKII_family / paras.v_tot;
end