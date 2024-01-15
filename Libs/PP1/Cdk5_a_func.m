function [Cdk5_a] = Cdk5_a_func(c_Ca, paras)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% v_Ca = c_Ca / (c_Ca + 0.5);
k0 = 0.1;
v_Ca = 0.1;
Cdk5_a = paras.Cdk5_tot / (1 + v_Ca/k0);
% v_Ca = 1;
% k_Ca = 0 / (1 + v_Ca);
% Cdk5_a = paras.Cdk5_tot * k_Ca;
end

