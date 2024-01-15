function [n_NMDAR_new] = Update_NMDAR(n_NMDAR, i_NMDAR)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% n_NMDAR / n_NMDAR_new -- number of NMDARs in different states: 4*1 array
% i_NMDAR -- index of the reaction for NMDARs: integer from 1 to 5

Update_lst = [1, -1, 0, 0;
              0, 0, 1, -1;
              1, 0, -1, 0;
              0, 1, 0, -1;
              0, -1, 0, 1]';

n_NMDAR_new = n_NMDAR + Update_lst(:, i_NMDAR);

end

