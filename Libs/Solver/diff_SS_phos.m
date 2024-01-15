function [d, S, eta_S_new] = diff_SS_phos(c_Ms, c_MPs, paras, c_Ca, c_E0, c_EP0)
    n = paras.n;

    [p_M, p_M_ATP, p_M_CaM, p_MA_CaM, p_M_EP, ...
     p_MP, p_MP_ADP, p_MP_CaM, p_MPA_CaM, p_MP_E, c_CaM4, c_E, c_EP] ...
     = compute_Prob(c_Ms, c_MPs, c_Ca, c_E0, c_EP0, paras.c_CaM, paras);

    gamma_c_E0 = (c_E0 + paras.c_E0_base) / c_E0;
    p_MP_E = p_MP_E * gamma_c_E0;

    pprate_nop_scale = (p_M_CaM+p_MA_CaM) * p_MA_CaM * paras.alpha_AO;
    pprate_nop_rev_scale = (p_M_CaM+p_MA_CaM) * p_MPA_CaM * paras.beta_AO;
    pprate_p_scale = (p_MP+p_MP_ADP+p_MP_CaM+p_MPA_CaM+p_MP_E) * p_MA_CaM * paras.alpha_AP;
    pprate_p_rev_scale = (p_MP+p_MP_ADP+p_MP_CaM+p_MPA_CaM+p_MP_E) * p_MPA_CaM * paras.beta_AP;
    sprate_scale = p_M_ATP * paras.alpha_A + p_MA_CaM * paras.alpha_AC + p_M_EP * paras.beta_E;
    sprate_rev_scale = p_MP_ADP * paras.beta_A + p_MPA_CaM * paras.beta_AC + p_MP_E * paras.alpha_E;
    
    S = pprate_nop_scale * paras.pprate_nop + pprate_nop_rev_scale * paras.pprate_nop_rev ...
      + pprate_p_scale * paras.pprate_p + pprate_p_rev_scale * paras.pprate_p_rev ...
      + sprate_scale * paras.sprate + sprate_rev_scale * paras.sprate_rev;
    
    eta_S = null(S); eta_S = eta_S(:, 1);
    eta_S = eta_S / sum(eta_S);
    
    eta_S_new = eta_S;
    eta_S_new(1) = (n*paras.c0 / (c_Ms+c_MPs) - 1) * paras.mu / (paras.nu * (p_M+p_M_ATP)^n);

    d(1) = (c_Ms + c_MPs) * dot(eta_S_new, paras.M_lst) - paras.n * c_Ms;
    d(2) = norm(S * eta_S_new);
end