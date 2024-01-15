function [n_CaV13_new] = Update_CaV13_channel(n_CaV13, i_CaV13)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% n_CaV13 / n_CaV13_new -- number of CaV13 channels in different states: 7*1 array
% i_CaV13 -- index of the reaction for CaV13 channels: integer from 1 to 20

Update_lst = [-1, 1, 0, 0, 0, 0, 0;
              1, -1, 0, 0, 0, 0, 0;
              1, 0, -1, 0, 0, 0, 0;
              -1, 0, 1, 0, 0, 0, 0;
              -1, 0, 0, 1, 0, 0, 0;
              1, 0, 0, -1, 0, 0, 0;
              0, 1, -1, 0, 0, 0, 0;
              0, -1, 1, 0, 0, 0, 0;
              0, 0, -1, 1, 0, 0, 0;
              0, 0, 1, -1, 0, 0, 0;
              0, -1, 0, 0, 1, 0, 0;
              0, 1, 0, 0, -1, 0, 0;
              0, 0, 1, 0, 0, -1, 0;
              0, 0, -1, 0, 0, 1, 0;
              0, 0, 0, -1, 0, 0, 1;
              0, 0, 0, 1, 0, 0, -1;
              0, 0, 0, 0, -1, 1, 0;
              0, 0, 0, 0, 1, -1, 0;
              0, 0, 0, 0, 0, 1, -1;
              0, 0, 0, 0, 0, -1, 1]';

n_CaV13_new = n_CaV13 + Update_lst(:, i_CaV13);

end

