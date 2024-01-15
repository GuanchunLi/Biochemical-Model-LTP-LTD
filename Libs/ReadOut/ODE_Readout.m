function [dy_lst, c_lst, pM_lst, pMP_lst] = ODE_Readout(Reaction_ode, t_lst, y_lst)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[m, n] = size(y_lst);
dy_lst = zeros(m, n);
c_lst = zeros(m, 2);
pM_lst = zeros(m, 5);
pMP_lst = zeros(m, 5);
for k = 1:m
    [dy_temp, c_temp, pM_temp, pMP_temp] = Reaction_ode(t_lst(k), y_lst(k, :)');
    dy_lst(k, :) = dy_temp';
    c_lst(k, :) = c_temp;
    pM_lst(k, :) = pM_temp;
    pMP_lst(k, :) = pMP_temp;
end

