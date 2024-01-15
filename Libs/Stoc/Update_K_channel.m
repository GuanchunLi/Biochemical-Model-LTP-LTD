function [n_K_new] = Update_K_channel(n_K, i_K)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% n_K / n_K_new -- number of potassium channels in different states: 5*1 array
% i_K -- index of the reaction for potassium channels: integer from 1 to 8

Update_lst = [-1, 1, 0, 0, 0;
              1, -1, 0, 0, 0;
              0, -1, 1, 0, 0;
              0, 1, -1, 0, 0;
              0, 0, -1, 1, 0;
              0, 0, 1, -1, 0;
              0, 0, 0, -1, 1;
              0, 0, 0, 1, -1]';

n_K_new = n_K + Update_lst(:, i_K);

end

