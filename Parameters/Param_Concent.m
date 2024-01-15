function [paras] = Param_Concent(paras)
%UNTITLED7 Summary of this function goes here
% Concentration of ATP family (out of equalibrium)
% paras.c_ATP = 1.0; paras.c_ADP = 0.1; paras.c_P = 0.1; 
paras.c_ATP = 100; paras.c_ADP = 0.1; paras.c_P = 0.1;

% paras.c_ADP = 0.1; paras.c_P = 1; % TEST CASE
% Concentration of all CaMi
paras.c_CaM = 30;
% paras.c_CaM = 24; % non-neigh-phos

% Concentration of CaMKII enzymes (6-unit)
paras.c0 = 12;

% Concentration of all enzyme molecular
% paras.c_E0 = 0.8;
paras.c_E_tot = 100;
% paras.c_E_tot = 1;

% NEW ADDED CODE
paras.c_MPs_base = 0;
paras.c_MP_base = 6e-8;
paras.c_E0_base = 1; % Should be 1
end

