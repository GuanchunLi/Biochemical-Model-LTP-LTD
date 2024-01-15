function [n_Na_new] = Update_Na_channel(n_Na, i_Na)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% n_Na / n_Na_new -- number of sodium channels in different states: 5*1 array
% i_Na -- index of the reaction for sodium channels: integer from 1 to 20

Update_lst = [-1, 1, 0, 0, 0, 0, 0, 0;
              1, -1, 0, 0, 0, 0, 0, 0;
              0, -1, 1, 0, 0, 0, 0, 0;
              0, 1, -1, 0, 0, 0, 0, 0;
              0, 0, -1, 1, 0, 0, 0, 0;
              0, 0, 1, -1, 0, 0, 0, 0;
              0, 0, 0, 0, -1, 1, 0, 0;
              0, 0, 0, 0, 1, -1, 0, 0;
              0, 0, 0, 0, 0, -1 ,1, 0;
              0, 0, 0, 0, 0, 1, -1, 0;
              0, 0, 0, 0, 0, 0, -1, 1;
              0, 0, 0, 0, 0, 0, 1, -1;
              -1, 0, 0, 0, 1, 0, 0, 0;
              1, 0, 0, 0, -1, 0, 0, 0;
              0, -1, 0, 0, 0, 1, 0, 0;
              0, 1, 0, 0, 0, -1, 0, 0;
              0, 0, -1, 0, 0, 0, 1, 0;
              0, 0, 1, 0, 0, 0, -1, 0;
              0, 0, 0, -1, 0, 0, 0, 1;
              0, 0, 0, 1, 0, 0, 0, -1]';

n_Na_new = n_Na + Update_lst(:, i_Na);

end

