function [c_Ms, c_MPs] = State2Con(state, paras)
%UNTITLED 此处提供此函数的摘要
%   此处提供详细说明
state_CaMKII = state(:, 1:end-1);
state_aCaMKII = state_CaMKII(:, state_CaMKII(1, :) == 1);
c_Ms = sum(sum(1 - state_aCaMKII(2:(paras.n+1), :))) / paras.coef_V;
c_MPs = sum(sum(state_aCaMKII(2:(paras.n+1), :))) / paras.coef_V;
end