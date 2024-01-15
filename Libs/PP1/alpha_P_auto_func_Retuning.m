function [alpha_P_auto_Ca] = alpha_P_auto_func_Retuning(c_Ca, paras)
% alpha_P_auto_Ca = paras.alpha_P_auto * (1 + 2 * (c_Ca > 1.4));
K = 0.24; n = 16;
% alpha_P_auto_Ca = paras.alpha_P_auto * (2000 * (c_Ca.^n) ./ (c_Ca.^n + K.^n));
alpha_P_auto_Ca = paras.alpha_P_auto * 0.01; % min(20, 2000 * (c_Ca.^n) ./ (c_Ca.^n + K.^n));
% alpha_P_auto_Ca = paras.alpha_P_auto * min(20, 3e-10*exp(150*c_Ca));
end