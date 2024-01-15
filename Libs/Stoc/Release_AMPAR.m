function [n_AMPAR_new] = Release_AMPAR(n_AMPAR)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% n_AMPAR / n_AMPAR_new -- number of AMPARs in different states: 4*1 array

n_AMPAR_new = n_AMPAR;
n_AMPAR_new(1) = 0;
n_AMPAR_new(2) = n_AMPAR(1) + n_AMPAR(2);
n_AMPAR_new(3) = 0;
n_AMPAR_new(4) = n_AMPAR(3) + n_AMPAR(4);

end

