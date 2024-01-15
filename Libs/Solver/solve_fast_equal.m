function X = solve_fast_equal(c_CaM, c_E0, c_EP0, c_Ms, c_MPs, K_Ca, paras)
func_equal = @(c_CaM4) diff_fast_equal_C4(c_CaM4,...
            c_CaM, c_E0, c_EP0, c_Ms, c_MPs, K_Ca, paras);

% options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective', 'SpecifyObjectiveGradient',true,...
%    'OptimalityTolerance', 1e-8, 'FunctionTolerance', 1e-8, 'Display', 'off'); 
% options = optimoptions('lsqnonlin', 'MaxIterations', 10000, 'StepTolerance', 1e-10, ...
%      'OptimalityTolerance', 1e-10, 'FunctionTolerance', 1e-10, 'Display', 'off'); 
% options = optimoptions('lsqnonlin', 'Display', 'off');

c_CaM4 = fzero(func_equal, [0, c_CaM]);

X = zeros(3, 1);
X(1) = c_CaM4;

k_ATP = 1 + paras.c_ATP / paras.K_MCaM_MACaM;
k_ADP = 1 + paras.c_ADP / paras.K_MPCaM_MPACaM;

d_c_M_func = @(c_M) c_Ms - c_M .* (1 + paras.c_ATP/paras.K_M_ATP + ...
    c_EP0./(paras.K_M_EP+c_M) + c_CaM4*k_ATP./paras.K_M_CaM);
c_M = fzero(d_c_M_func, [0, c_Ms]);

d_c_MP_func = @(c_MP) c_MPs - c_MP .* (1 + paras.c_ADP/paras.K_MP_ADP + ...
    c_E0./(paras.K_MP_E+c_MP) + c_CaM4*k_ADP./paras.K_MP_CaM);
c_MP = fzero(d_c_MP_func, [0, c_MPs]);

X(2) = c_M; X(3) = c_MP;

end

function d = diff_fast_equal_C4(c_CaM4,...
            c_CaM, c_E0, c_EP0, c_Ms, c_MPs, K_Ca, paras)

k_ATP = 1 + paras.c_ATP / paras.K_MCaM_MACaM;
k_ADP = 1 + paras.c_ADP / paras.K_MPCaM_MPACaM;

d_c_M_func = @(c_M) c_Ms - c_M .* (1 + paras.c_ATP/paras.K_M_ATP + ...
    c_EP0./(paras.K_M_EP+c_M) + c_CaM4*k_ATP./paras.K_M_CaM);
c_M = fzero(d_c_M_func, [0, c_Ms]);

d_c_MP_func = @(c_MP) c_MPs - c_MP .* (1 + paras.c_ADP/paras.K_MP_ADP + ...
    c_E0./(paras.K_MP_E+c_MP) + c_CaM4*k_ADP./paras.K_MP_CaM);

c_MP = fzero(d_c_MP_func, [0, c_MPs]);

d = c_CaM - c_CaM4 * (K_Ca + c_M*k_ATP./paras.K_M_CaM + c_MP*k_ADP./paras.K_MP_CaM);
end
