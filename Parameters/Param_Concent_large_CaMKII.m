function [paras] = Param_Concent_large_CaMKII(paras)
%UNTITLED7 Summary of this function goes here
% Concentration of ATP family (out of equalibrium)
% paras.c_ATP = 1.0; paras.c_ADP = 0.1; paras.c_P = 0.1; 
paras.c_ATP = 1000; paras.c_ADP = 0.0; paras.c_P = 0.0;

% Concentration of all CaMi
paras.c_CaM = 30;
% paras.c_CaM = 24; % non-neigh-phos

% Concentration of all enzyme molecular
paras.c_E0 = 0.1;
end

