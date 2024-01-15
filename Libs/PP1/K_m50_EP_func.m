function [K_m50_EP] = K_m50_EP_func(c_Ca, paras)
% alpha_P_Ca = paras.alpha_P + 0.18 * (c_Ca > 1.4); % 0.5 * (c_Ca.^3) ./ (c_Ca.^3 + 1.4.^3);
% n = 6; % 4 -- original value
% K_Ca = 1 + 8e-5 / (8e-12 + 0*max(c_Ca-0.1) .^ n);
% c_Ca = min(c_Ca, 0.18);
% K_Ca = 1e7 * exp(-128*(c_Ca-0.1)) + 1;
K_Ca = 1 + 1e6*exp(-92*(c_Ca-0.1));
K_m50_EP = paras.K_m50_EP * K_Ca;
end