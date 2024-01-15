function [alpha_P_Ca] = alpha_P_func_Retuning(c_Ca, paras)
K = 0.1; n = 4; % 4 -- original value
% alpha_P_Ca = paras.alpha_P * (c_Ca - 0.1) .^ n / ((c_Ca - 0.1).^n + K.^n);
% v_Ca = max(c_Ca-0.1, 0).^2 ./ (max(c_Ca-0.1, 0).^2 + 0.09.^2) + 5e-3; Originall
v_Ca = max(c_Ca-0.1, 0).^2 ./ (max(c_Ca-0.1, 0).^2 + 0.09.^2) + 5e-3;
alpha_P_Ca = paras.alpha_P * v_Ca;
end