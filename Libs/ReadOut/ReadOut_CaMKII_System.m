function [K_MP_E_1, K_MP_E_2, p_MP] ...
    = ReadOut_CaMKII_System(t, y_CaMKII, c_Ca, c_E0, c_EP0, c_CaM, paras)
% READOUT_CAMKII_SYSTEM
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
    gamma_c_MPs = (c_MPs + paras.c_MPs_base) / c_MPs;
    gamma_c_E0 = (c_E0 + paras.c_E0_base) / c_E0;
    
    dE_CaM = p_M_EP * paras.beta_E * c_Ms - p_MP_E * paras.alpha_E * c_MPs * gamma_c_MPs;
    
    % Compute the avaliable active CaM
    K_Ca = equal_constant_Ca(c_Ca, paras);
    c_CaM4_tot = c_CaM - (K_Ca - 1) * c_CaM4;
    % c_CaM4_tot = c_CaM4; % TEST CASE
    
    K_m_CaM4 = 1 + paras.c_ADP/paras.K_MP_ADP + ...
        c_CaM4/paras.K_MP_CaM * (1 + paras.c_ADP/paras.K_MPCaM_MPACaM);
    
   K_MP_E_1 = paras.K_MP_E/p_MP;
   K_MP_E_2 = K_m_CaM4 * (paras.K_MP_E + p_MP*c_MPs);
end