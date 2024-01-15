function [n_NMDAR_new] = Release_NMDAR(n_NMDAR)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% n_NMDAR / n_NMDAR_new -- number of AMPARs in different states: 4*1 array

n_NMDAR_new = n_NMDAR;
n_NMDAR_new(1) = 0;
n_NMDAR_new(2) = n_NMDAR(1) + n_NMDAR(2);
n_NMDAR_new(3) = 0;
n_NMDAR_new(4) = n_NMDAR(3) + n_NMDAR(4);

end

