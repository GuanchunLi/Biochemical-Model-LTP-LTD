function [v_Ca] = Cdk5_Ca_func(c_Ca, paras)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% mod 4 v_Ca = c_Ca.^6 ./ (c_Ca.^6 + 0.5.^6);
% v_Ca = 1;
v_Ca = 2e5 * (c_Ca-0.1).^6 ./ ((c_Ca-0.1).^6 + 0.4.^6);
v_Ca = v_Ca / paras.v0;
end

