function [X, flag] = solve_SS_phos(c_Ca, c_E0, c_EP0, paras)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
func_equal = @(X) diff_SS_phos(X(1), X(2), paras, c_Ca, c_E0, c_EP0);
c_Ms_range = [0, paras.n*paras.c0]; c_MPs_range = [0, paras.n*paras.c0];
[X, flag] = solve_range(func_equal, c_Ms_range, c_MPs_range, paras.n*paras.c0);
while max(flag) < 0
    [X, flag] = solve_range(func_equal, c_Ms_range, c_MPs_range, paras.n*paras.c0);
end
end

