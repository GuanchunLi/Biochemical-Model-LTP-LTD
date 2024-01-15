function [n_CaV12_new] = Update_CaV12_channel(n_CaV12, i_CaV12)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% n_CaV12 / n_CaV12_new -- number of CaV13 channels in different states: 7*1 array
% i_CaV12 -- index of the reaction for CaV13 channels: integer from 1 to 20

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

n_CaV12_new = n_CaV12 + Update_lst(:, i_CaV12);

end

