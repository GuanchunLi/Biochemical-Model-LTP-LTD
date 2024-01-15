function [n_AMPAR_new] = Update_AMPAR(n_AMPAR, i_AMPAR)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% n_AMPAR / n_AMPAR_new -- number of AMPARs in different states: 4*1 array
% i_AMPAR -- index of the reaction for AMPARs: integer from 1 to 5

Update_lst = [1, -1, 0, 0;
              0, 0, 1, -1;
              1, 0, -1, 0;
              0, 1, 0, -1;
              0, -1, 0, 1]';

n_AMPAR_new = n_AMPAR + Update_lst(:, i_AMPAR);

end

