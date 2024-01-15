function [paras] = Param_CaM_Linse(paras)
%PARAM_CAM Set parameters for the binding of Ca2+ to CaM

% % Linse(1991), [KCL] Low
% paras.K_1 = 0.1; % 0.06~0.16
% paras.K_2 = 0.251; % 0.16~0.04
% paras.K_3 = 0.316; % 0.2~0.5
% paras.K_4 = 0.4; % 0.25~0.63

% Linse(1991), [KCL] = 25mM
paras.K_1 = 1.6; % 1.6, 1~2.5
paras.K_2 = 0.1; % 0.1, 0.063~0.16
paras.K_3 = 4; % 2.5~6.31
paras.K_4 = 2; % 1.25~3.16

% % Linse(1991), [KCL] = 50mM
% paras.K_1 = 4; % 2.5~6.3
% paras.K_2 = 0.2; % 0.125~0.32
% paras.K_3 = 25; % 16~40
% paras.K_4 = 3; % 2~5

% % Linse(1991), [KCL] = 150mM
% paras.K_1 = 20; % 12.6~31.6
% paras.K_2 = 0.63; % 0.4~1
% paras.K_3 = 100; % 50~200
% paras.K_4 = 5; % 2.5~10

end

