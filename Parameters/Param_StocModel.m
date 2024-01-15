function [paras] = Param_StocModel(paras, vol)
%PARAM_RATEMODEL Aux variables for the stochastic model
paras.reaction_lst = [1; (2:(paras.n+1))'];
paras.coef_V = 6.022 * 10^2 * vol;
end

