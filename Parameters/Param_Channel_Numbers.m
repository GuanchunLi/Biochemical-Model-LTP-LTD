function [paras] = Param_Channel_Numbers(paras)
%UNTITLED3 此处提供此函数的摘要
%   此处提供详细说明
paras.n_AMPAR = 120 * paras.coef_V_channel;
paras.n_NMDAR = 7 * paras.coef_V_channel;
paras.n_CaV12 = 4 * 1 * paras.coef_V_channel;
paras.n_CaV13 = 1 * 1 * paras.coef_V_channel;
% % An area of 5 muM^2
% paras.n_Na = 42; % 42 ; % 320 pS/muM^2, single channel 17pS
% paras.n_K = 185; % 185; % 480 pS/muM^2, single channel 7pS
end