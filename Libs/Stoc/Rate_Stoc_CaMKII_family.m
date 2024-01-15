function [R_mat_all] = Rate_Stoc_CaMKII_family(t, y, c_Ca, paras)
% The reaction scheme for the CaMKII subunits with neighboring
% phosphorylation (Stochastic version)
% INPUTS:
%        t -- s Current time 
%        y -- a (1+6) * (n_CaMKII+1) matrix
%        c_Ca -- muM Concentration of Calcium at current time     
%        paras -- Group of parameters determing the fast equalibrium
% OUTPUTS:
%        R_mat_all -- a (1+5*6) * (1+n_CaMKII) matrix
% Dependency: 
%        compute_Prob_special.m

    % c_Ca = cal_source(t);
    n = paras.n;
    c_E0 = y(1, end) / paras.coef_V; c_EP0 = y(2, end) / paras.coef_V;
    y_CaMKII = y(:, 1:end-1);
    y_aCaMKII = y_CaMKII(:, y_CaMKII(1, :) == 1);
    c_Ms = sum(sum(1 - y_aCaMKII(2:(n+1), :))) / paras.coef_V;
    c_MPs = sum(sum(y_aCaMKII(2:(n+1), :))) / paras.coef_V;
    
    [p_M, p_M_ATP, p_M_CaM, p_MA_CaM, p_M_EP, ...
     p_MP, p_MP_ADP, p_MP_CaM, p_MPA_CaM, p_MP_E, c_CaM4, c_E, c_EP] ...
     = compute_Prob(c_Ms, c_MPs, c_Ca, c_E0, c_EP0, paras.c_CaM, paras);
    
    % Deal with the case when c_E0 = 0
    if c_E0 > 0
        gamma_c_E0 = (c_E0 + paras.c_E0_base) / c_E0;
    else
        gamma_c_E0 = 1;
    end
    p_MP_E_base = paras.c_E0_base * p_MP / (p_MP * c_MPs + paras.K_MP_E) * 1;

    c_MP_E_base_dephos = c_E * paras.c_MP_base / paras.K_MP_E;

    % Reaction rate for each reaction
    [nx, ny] = size(y); n_CaMKII = ny-1;
    R_mat_all = zeros(1+6*n, n_CaMKII+1);
    % Dynamics of PP
    % Auto-dephosphorylation
    R_mat_all(1, end) = (paras.alpha_P_auto * c_EP / (paras.K_m50_EP + c_EP) * c_E ...
                         - paras.alpha_E * c_MP_E_base_dephos) * paras.coef_V;
    % Calcium-induced dephosphorylation
    alpha_P_Ca = alpha_P_func_Retuning(c_Ca, paras);
    R_mat_all(2, end) = (alpha_P_Ca * c_EP + 0 * paras.d_E_base) * paras.coef_V;

    % Dynamics of CaMKII
    % Opening and closing of CaMKII
    mu = paras.mu * 5 * (c_Ca^2) / (c_Ca^2+0.2^2);
    R_mat_all(1, 1:end-1) = mu * (1 - y_CaMKII(1, :)) + paras.nu * ((sum(y_CaMKII(2:end, :), 1) < 0.5) .* y_CaMKII(1, :)) * (p_M+p_M_ATP)^n;
    % Self phosphorylation without CaM
    R_mat_all(2:(n+1), 1:end-1) = y_CaMKII(1, :) .* p_MP_ADP * paras.beta_A .* y_CaMKII(2:(n+1), :)...
                          + y_CaMKII(1, :) .* p_M_ATP * paras.alpha_A .* (1 - y_CaMKII(2:(n+1), :));
    % Self phosphorylation with CaM
    R_mat_all((n+2):(2*n+1), 1:end-1) = y_CaMKII(1, :) .* p_MPA_CaM * paras.beta_AC .* y_CaMKII(2:(n+1), :) ...
                                + y_CaMKII(1, :) .* p_MA_CaM * paras.alpha_AC .* (1 - y_CaMKII(2:(n+1), :));
    % Neighbouring phosphorylation (catalyst with CaM)
    R_mat_all((2*n+2):(3*n+1), 1:end-1) = y_CaMKII(1, :) .* p_MA_CaM * (p_MA_CaM + p_M_CaM) * paras.alpha_AO .* (1 - y_CaMKII(2:(n+1), :)) .* (1 - [y_CaMKII(3:(n+1), :); y_CaMKII(2, :)])... 
                                  + y_CaMKII(1, :) .* p_MPA_CaM * (p_MA_CaM + p_M_CaM) * paras.beta_AO .* y_CaMKII(2:(n+1), :) .* (1 - [y_CaMKII(3:(n+1), :); y_CaMKII(2, :)]);
    % Neighbouring phosphorylation (catalyst phosphorylated)
    R_mat_all((3*n+2):(4*n+1), 1:end-1) = y_CaMKII(1, :) .* (p_MP+p_MP_ADP+p_MP_CaM+p_MPA_CaM) * p_MA_CaM * paras.alpha_AP .* (1 - y_CaMKII(2:(n+1), :)) .* [y_CaMKII(3:(n+1), :); y_CaMKII(2, :)]... 
                                  + y_CaMKII(1, :) .* (p_MP+p_MP_ADP+p_MP_CaM+p_MPA_CaM) * p_MPA_CaM * paras.beta_AP .* y_CaMKII(2:(n+1), :) .* [y_CaMKII(3:(n+1), :); y_CaMKII(2, :)];
    if c_E0 > 0
    % Dephosphoryolation with Baseline PP
    R_mat_all((4*n+2):(5*n+1), 1:end-1) = y_CaMKII(1, :) .* p_MP_E * (gamma_c_E0 - 1) * paras.alpha_E .* y_CaMKII(2:(n+1), :)... 
                                  + y_CaMKII(1, :) .* p_M_EP * (gamma_c_E0 - 1) * paras.beta_E .* (1 - y_CaMKII(2:(n+1), :));
    else
        % Dephosphoryolation with Baseline PP
    R_mat_all((4*n+2):(5*n+1), 1:end-1) = y_CaMKII(1, :) .* p_MP_E_base * paras.alpha_E .* y_CaMKII(2:(n+1), :)... 
                                  + y_CaMKII(1, :) .* p_M_EP * (gamma_c_E0 - 1) * paras.beta_E .* (1 - y_CaMKII(2:(n+1), :));
    end
    % Dephosphoryolation with active PP
    R_mat_all((5*n+2):(6*n+1), 1:end-1) = y_CaMKII(1, :) .* p_MP_E * paras.alpha_E .* y_CaMKII(2:(n+1), :)... 
                                  + y_CaMKII(1, :) .* p_M_EP * paras.beta_E .* (1 - y_CaMKII(2:(n+1), :));
        
    % t_mat = (-log(rand(size(R_mat)))) ./ R_mat;
%     t_mat = (-log(rand(size(R_mat_all)))) ./ R_mat_all;
%     [tk, ik] = min(abs(t_mat), [], 1);
%     [tj, jk] = min(tk, [], 2);
%     idx = [ik(jk), jk];
%     t_new = t + t_mat(ik(jk), jk);
    R_mat_all(R_mat_all < 0) = 0;
    R_mat_all = R_mat_all / paras.v_tot;
end