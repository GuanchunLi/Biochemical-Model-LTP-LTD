function [AMPA_change] = AMPA_SS(r_MPs, r_PPs, paras)
%UNTITLED 此处提供此函数的摘要
% r_AMPA -- Ratio of Active AMPA Channels (linear with EPSP)
% 此处提供详细说明
k_AMPA_phos = paras.k_AMPA_phos0 + r_MPs * paras.d_AMPA_phos / (18*6);
k_AMPA_dephos = paras.k_AMPA_dephos0 + r_PPs * paras.d_AMPA_dephos / (100);
AMPA_change = (k_AMPA_phos ./ (k_AMPA_phos + k_AMPA_dephos)) ...
             / (paras.k_AMPA_phos0 / (paras.k_AMPA_phos0 + paras.k_AMPA_dephos0));
end