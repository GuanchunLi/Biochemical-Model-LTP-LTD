function [R_mat_all] = Rate_stoc_wi_transpp(t, y, c_Ca, c_E0, c_CaM, paras)
% The reaction scheme for the CaMKII subunits with neighboring
% phosphorylation (Stochastic version)
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
% Dependency: 
%        compute_Prob_special.m

    % c_Ca = cal_source(t);
    n = paras.n;
    y_o = y(:, y(1, :) == 1);
    c_Ms = sum(sum(1 - y_o(2:(n+1), :))) / paras.coef_V;
    c_MPs = sum(sum(y_o(2:(n+1), :))) / paras.coef_V;
    
    [p_M, p_M_ATP, p_M_CaM, p_MA_CaM, p_M_EP, ...
     p_MP, p_MP_ADP, p_MP_CaM, p_MPA_CaM, p_MP_E] ...
     = compute_Prob_special(c_Ms, c_MPs, c_Ca, c_E0, c_CaM, paras);
 
    % Reaction rate for each reaction
    [nx, ny] = size(y);
    R_mat = zeros(size(y));
    R_mat(1, :) = paras.mu * (1 - y(1, :)) + paras.nu * ((sum(y(2:end, :), 1) < 0.5) .* y(1, :)) * (p_M+p_M_ATP)^n;
    R_mat(2:(n+1), :) = y(1, :) .* p_MP_ADP * paras.beta_A .* y(2:(n+1), :)...
                      + y(1, :) .* p_M_ATP * paras.alpha_A .* (1 - y(2:(n+1), :));
    R_mat((n+2):(2*n+1), :) = y(1, :) .* p_MPA_CaM * paras.beta_AC .* y(2:(n+1), :) ...
                            + y(1, :) .* p_MA_CaM * paras.alpha_AC .* (1 - y(2:(n+1), :));
    R_mat((2*n+2):(3*n+1), :) = y(1, :) .* p_MP_E * paras.alpha_E .* y(2:(n+1), :)... 
                              + y(1, :) .* p_M_EP * paras.beta_E .* (1 - y(2:(n+1), :));
    % Neighbouring phosphorylation
    R_mat((3*n+2):(4*n+1), :) = y(1, :) .* p_MA_CaM * (p_MA_CaM + p_M_CaM) * paras.alpha_AO .* (1 - y(2:(n+1), :)) .* (1 - [y(3:(n+1), :); y(2, :)])... 
                              + y(1, :) .* p_MPA_CaM * (p_MA_CaM + p_M_CaM) * paras.beta_AO .* y(2:(n+1), :) .* (1 - [y(3:(n+1), :); y(2, :)]);
    R_mat((4*n+2):(5*n+1), :) = y(1, :) .* p_MA_CaM * paras.alpha_AP .* (1 - y(2:(n+1), :)) .* [y(3:(n+1), :); y(2, :)]... 
                              + y(1, :) .* p_MPA_CaM * paras.beta_AP .* y(2:(n+1), :) .* [y(3:(n+1), :); y(2, :)];
    R_mat_all = zeros(n+1, ny);
    R_mat_all(1, :) = R_mat(1, :);
    R_mat_all(2:(n+1), :) = R_mat(2:(n+1), :) + R_mat((n+2):(2*n+1), :) + R_mat((2*n+2):(3*n+1), :) + R_mat((3*n+2):(4*n+1), :) + R_mat((4*n+2):(5*n+1), :);
    % t_mat = (-log(rand(size(R_mat)))) ./ R_mat;
    t_mat = (-log(rand(size(R_mat_all)))) ./ R_mat_all;
    [tk, ik] = min(abs(t_mat), [], 1);
    [tj, jk] = min(tk, [], 2);
    idx = [ik(jk), jk];
    t_new = t + t_mat(ik(jk), jk);
    
%     fprintf(['t = ', num2str(t), ' c_Ms = ', num2str(c_Ms), ' c_MPs = ', num2str(c_MPs),...
%              ' c_MP = ', num2str(c_MPs*p_MP), '\n']);
end