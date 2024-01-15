function [dy_CaMKII, c_lst, pM_lst, pMP_lst, p_E, p_EP, dE_CaM, c_CaM4_tot] ...
    = Reaction_ode_wi_transpp_RE(t, y_CaMKII, c_Ca, c_E0, c_EP0, c_CaM, paras)
% REACTION_ODE_WI_TRANSPP
% The reaction scheme for the CaMKII subunits with neighboring
% phosphorylation
% INPUTS:
%        t -- s Current time 
%        y_CaMKII -- muM Concentration of CaMKII in each state
%        c_Ca -- muM Concentration of Calcium at current time
%        c_E0 -- muM Concentration of avaliable phosphatase
%        c_CaM -- muM Concentration of avaliable calmodulin        
%        paras -- Group of parameters determing the fast equalibrium
% OUTPUTS:
%        dy_CaMKII -- muM/s Change rate of CaMKII in each state
%        c_lst -- muM Concentration of unphosphorylated and phophorylated
%        CaMKII subunits
%        pM_lst -- 1 Probabilities of states for unphosphorylated subunits
%        pMP_lst -- 1 Probabilities of states for phosphorylated subunits
%        p_E -- 1 Probability of PP1 in the free state
%        c_CaM4 -- muM Concentration of available calmodulin with 4 Ca2+ ions
% Dependency: 
%        compute_Prob_special_new.m, equal_constant_Ca.m

    % c_Ca = cal_source(t);
    n = paras.n;
    y_open = y_CaMKII(2:end);
    c_Ms = dot(paras.M_lst, y_open);
    c_MPs = dot(paras.MP_lst, y_open);
    c_lst = [c_Ms, c_MPs];
    
%     % Baseline activity of PP1 added
%     c_E_base = 0.8;
%     c_E0 = c_E0 + c_E_base;
    
    [p_M, p_M_ATP, p_M_CaM, p_MA_CaM, p_M_EP, ...
     p_MP, p_MP_ADP, p_MP_CaM, p_MPA_CaM, p_MP_E, c_CaM4, c_E, c_EP] ...
     = compute_Prob(c_Ms, c_MPs, c_Ca, c_E0, c_EP0, c_CaM, paras);
     % = compute_Prob_special(c_Ms, c_MPs, c_Ca, c_E0, c_CaM, paras);
     % = compute_Prob(c_Ms, c_MPs, c_Ca, c_E0, c_CaM, paras);
    
    % Compute the probability of free PP1
    % p_E = (c_E - paras.c_E0_base) / c_E0;
    p_E = c_E  / c_E0;
    p_EP = c_EP / c_EP0;
    % Compute the reaction rate of PP1 due to dephosphorylating CaMKII
    
    % ADDED CODE
    % gamma_c_MPs = (c_MPs + paras.c_MPs_base) / c_MPs;
    c_MP = c_MPs*p_MP;
    gamma_c_MPs = 1 + paras.c_MP_base / c_MP;
    % gamma_c_MPs = 1;
    gamma_c_E0 = (c_E0 + paras.c_E0_base) / c_E0;

    dE_CaM = p_M_EP * paras.beta_E * c_Ms - p_MP_E * paras.alpha_E * c_MPs * gamma_c_MPs;
    
    % Compute the avaliable active CaM
    K_Ca = equal_constant_Ca(c_Ca, paras);
    c_CaM4_tot = c_CaM - (K_Ca - 1) * c_CaM4;
    % c_CaM4_tot = c_CaM4; % TEST CASE
    
    
    pM_lst = [p_M, p_M_ATP, p_M_CaM, p_MA_CaM, p_M_EP];
    pMP_lst = [p_MP, p_MP_ADP, p_MP_CaM, p_MPA_CaM, p_MP_E];
    
    % ADDED CODE
    p_MP_E = p_MP_E * gamma_c_E0;
    p_M_EP = p_M_EP * gamma_c_E0; % Modified for detailed balance
    
    % Reaction rate for each reaction
    pprate_nop_scale = (p_M_CaM+p_MA_CaM) * p_MA_CaM * paras.alpha_AO;
    pprate_nop_rev_scale = (p_M_CaM+p_MA_CaM) * p_MPA_CaM * paras.beta_AO;
    pprate_p_scale = (p_MP+p_MP_ADP+p_MP_CaM+p_MPA_CaM) * p_MA_CaM * paras.alpha_AP;
    pprate_p_rev_scale = (p_MP+p_MP_ADP+p_MP_CaM+p_MPA_CaM) * p_MPA_CaM * paras.beta_AP;
    sprate_scale = p_M_ATP * paras.alpha_A + p_MA_CaM * paras.alpha_AC + p_M_EP * paras.beta_E;
    sprate_rev_scale = p_MP_ADP * paras.beta_A + p_MPA_CaM * paras.beta_AC + p_MP_E * paras.alpha_E;
    
    % Calcium dependent opening
    mu = paras.mu * 5 * (c_Ca^2) / (c_Ca^2+0.2^2);
    dy_CaMKII = 0 * y_CaMKII;

    dy_CaMKII(2:end) = pprate_nop_scale * (paras.pprate_nop * y_open) + pprate_nop_rev_scale * (paras.pprate_nop_rev * y_open) +...
        + pprate_p_scale * (paras.pprate_p * y_open) + pprate_p_rev_scale * (paras.pprate_p_rev * y_open) +...
        + sprate_scale * (paras.sprate * y_open) + sprate_rev_scale * (paras.sprate_rev * y_open);
    dy_CaMKII(1) = - mu * y_CaMKII(1) + paras.nu * (p_M+p_M_ATP)^n * y_open(1);
    dy_CaMKII(2) = dy_CaMKII(2) + mu * y_CaMKII(1) - paras.nu * (p_M+p_M_ATP)^n * y_open(1);
%    fprintf(['t = ', num2str(t), ' c_Ms = ', num2str(c_Ms), ' c_MPs = ', num2str(c_MPs),...
%             ' c_MP = ', num2str(c_MPs*p_MP), ' p_Ms = ', num2str(sum(pM_lst)), ' p_MA_CaM = ', num2str(p_MA_CaM),'\n']);
end