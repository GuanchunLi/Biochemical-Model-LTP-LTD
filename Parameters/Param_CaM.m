function [paras] = Param_CaM(paras)
%PARAM_CAM Set parameters for the binding of Ca2+ to CaM
% r1 = 0.2; r2 = 0.1; v = 0.7;
% K_1_C = 2 * v; % 7.3; % 7.3 - 12
% K_2_C = K_1_C * r1 * v; % 0.4 - 1.7
% K_1_N = 5 * v; % 4;   % 15 - 40
% K_2_N = K_1_N * r2 * v; % 3.3; % 3.3- 9

% K_1_C = 10; % 7.3 - 12
% K_2_C = 1; % 0.4 - 1.7
% K_1_N = 20;   % 15 - 40
% K_2_N = 5; % 3.3- 9

% paras.K_4 = K_2_N + K_2_C;
% paras.K_3 = (K_1_C*K_2_C + K_2_N*K_2_C + K_1_N*K_2_N) / paras.K_4;
% paras.K_2 = (K_2_N*K_1_C*K_2_C + K_1_N*K_2_N*K_2_C) / paras.K_3;
% paras.K_1 = K_1_N*K_2_N*K_1_C*K_2_C / paras.K_2;

paras.K_1 = 0.1;
paras.K_2 = 0.02;
paras.K_3 = 5;
paras.K_4 = 5;

% % a4 = 0.025; a3 = 0.5; a2 = 40; a1 = 40; 
% a4 = 1; a3 = 0.2; a2 = 10; a1 = 5; 
% paras.K_1 = a4 / a3;
% paras.K_2 = a3 / a2;
% paras.K_3 = a2 / a1;
% paras.K_4 = a1;
end

