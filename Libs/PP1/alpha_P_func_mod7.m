function [alpha_P_Ca] = alpha_P_func_mod7(c_Ca, paras)
% K = 0.1; n = 8; % 4 -- original value
% P_Ca = (c_Ca-0.1).^n / ((c_Ca-0.1).^n + K.^n);
% c_Ca = min(c_Ca, 0.18);
P_Ca = 1e-6*exp(147.5*(c_Ca-0.1));
P_Ca = min(P_Ca, 0.5);
alpha_P_Ca = paras.alpha_P * P_Ca;
end