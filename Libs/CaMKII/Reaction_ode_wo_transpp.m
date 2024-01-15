function [dy_CaMKII, c_lst, pM_lst, pMP_lst, p_E] = Reaction_ode_wo_transpp(t, y_CaMKII, c_Ca, c_E0, c_CaM, paras)
% REACTION_ODE_WO_TRANSPP
% The reaction scheme for the CaMKII subunits without neighboring
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
% Dependency: 
%        compute_Prob_special.m

    % c_Ca = cal_source(t);
    n = paras.n;
    y_open = y_CaMKII(2:end);
    c_Ms = dot(paras.M_lst, y_open);
    c_MPs = dot(paras.MP_lst, y_open);
    c_lst = [c_Ms, c_MPs];
    
    [p_M, p_M_ATP, p_M_CaM, p_MA_CaM, p_M_EP, ...
     p_MP, p_MP_ADP, p_MP_CaM, p_MPA_CaM, p_MP_E] ...
     = compute_Prob_special(c_Ms, c_MPs, c_Ca, c_E0, c_CaM, paras);
    
    % Compute the probability of free PP1
    p_E = (c_E0 - p_MP_E * c_MPs) / c_E0;
    
    pM_lst = [p_M, p_M_ATP, p_M_CaM, p_MA_CaM, p_M_EP];
    pMP_lst = [p_MP, p_MP_ADP, p_MP_CaM, p_MPA_CaM, p_MP_E];
    
    % Compute the unit rate
    gamma_u = paras.alpha_A * p_M_ATP + paras.alpha_AC * p_MA_CaM + paras.beta_E * p_M_EP;
    delta_u = paras.beta_A * p_MP_ADP + paras.beta_AC * p_MPA_CaM + paras.alpha_E * p_MP_E;

    dy_CaMKII = 0 * y_CaMKII;
    % Rate for the closed state
    dy_CaMKII(1) = - paras.mu * y_CaMKII(1) + paras.nu * (p_M+p_M_ATP)^n * y_open(1);
    
    % Rate for the O_0 state
    dy_CaMKII(2) = paras.mu * y_CaMKII(1) - (paras.nu * (p_M+p_M_ATP)^n + n * gamma_u) * y_open(1) + delta_u * y_open(2);
    
    % Rate for the O_1 - O_n-1 states
    for k = 1:(n-1)
        dy_CaMKII(k+2) = (n-k+1)*gamma_u * y_open(k) - ((n-k)*gamma_u + k*delta_u)*y_open(k+1) + (k+1)*delta_u * y_open(k+2);
    end
    
    % Rate for the O_n state
    dy_CaMKII(n+2) = gamma_u * y_open(n) - n*delta_u * y_open(n+1);
    
    fprintf(['t = ', num2str(t), ' c_Ms = ', num2str(c_Ms), ' c_MPs = ', num2str(c_MPs),...
             ' c_MP = ', num2str(c_MPs*p_MP), ' p_Ms = ', num2str(p_M+p_M_ATP), ' p_MA_CaM = ', num2str(p_MA_CaM),'\n']);
end