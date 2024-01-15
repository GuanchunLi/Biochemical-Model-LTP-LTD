function [paras] = Param_FastEqual_Retuning(paras)
%PARAM_FASTEQUAL Set parameters for the fast equalibrium constant
% Detailed balance: K_M_ATP * K_MA_MACaM = K_M_CaM * K_MCaM_MACaM
paras.K_M_ATP = 1;
paras.K_M_CaM = 1/(10); %1/10
paras.K_MCaM_MACaM = 1;
paras.K_MA_MACaM = 1/(10); %1/10

% Detailed balance: K_MP_ADP * K_MPA_MPACaM = K_MP_CaM * K_MPCaM_MPACaM
paras.K_MP_ADP = 1;
paras.K_MP_CaM = 1/(10^3); %1/10^3
paras.K_MPCaM_MPACaM = 1;
paras.K_MPA_MPACaM = 1/(10^3); %1/10^3

% Fast equalibrium with the phosphotase
paras.K_M_EP = 10^3; %10^3
paras.K_MP_E = 1/(10^3) * 1; %1/10^3
% paras.K_EP = 10^4;
end

