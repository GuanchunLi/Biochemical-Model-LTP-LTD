function [alpha_P_Ca] = alpha_P_func(c_Ca, paras)
K = 0.1; n = 4; % 4 -- original value
alpha_P_Ca = paras.alpha_P * (c_Ca - 0.1) .^ n / ((c_Ca - 0.1).^n + K.^n);
% alpha_P_Ca = paras.alpha_P * 5e-4;
end