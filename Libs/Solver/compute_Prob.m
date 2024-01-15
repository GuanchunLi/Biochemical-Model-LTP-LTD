function [p_M, p_M_ATP, p_M_CaM, p_MA_CaM, p_M_EP, ...
          p_MP, p_MP_ADP, p_MP_CaM, p_MPA_CaM, p_MP_E, ...
          c_CaM4, c_E, c_EP] ...
          = compute_Prob(c_Ms, c_MPs, c_Ca, c_E0, c_EP0, c_CaM, paras)
% COMPUTE_PROB
% Compute the probability of each state of CaMKII subunit 
% in the fast-equalibrium case 
% Inputs:
%        c_Ms -- muM Concentration of unphosphorylated subunits 
%        c_MPs -- muM Concentration of unphosphorylated subunits
%        c_Ca -- muM Concentration of Calcium
%        c_E0 -- muM Concentration of avaliable phosphatase
%        c_CaM -- muM Concentration of avaliable calmodulin
%        paras -- Group of parameters determing the fast equalibrium
% Outputs:
%        p_Mx -- 1 Probability of an unphosphorylated subunit in state Mx
%        p_MPx -- 1 Probability of a phosphorylated subunit in state MPx
%        c_CaM4 -- muM Concentration of Calmodulin with 4 Ca2+ ions
%        c_E -- muM Concentration of free phosphatase
% Dependency:
%        equal_constant_Ca.m,  solve_fast_equal.m


% Initialization
K_Ca = equal_constant_Ca(c_Ca, paras);

% Solve the system
X = solve_fast_equal(c_CaM, c_E0, c_EP0, c_Ms, c_MPs, K_Ca, paras);

c_CaM4 = X(1); c_M = X(2); c_MP = X(3);

% Compute the concentrations of all other chemicals
k_ATP = 1 + paras.c_ATP / paras.K_MCaM_MACaM;
k_ADP = 1 + paras.c_ADP / paras.K_MPCaM_MPACaM;
c_E = c_E0 / (1 + c_MP/paras.K_MP_E);
c_EP = c_EP0 / (1 + c_M/paras.K_M_EP);

% Throw error when the solution is non-reasonable
if (c_M < 0 || c_M > c_Ms || c_MP < 0 || c_MP > c_MPs || c_CaM4 < 0 || ...
    c_CaM4 > c_CaM || c_E < 0 || c_E > c_E0 || c_EP < 0 || c_EP > c_EP0)
    fprintf(['c_E = ', num2str(c_E), ' c_EP = ', num2str(c_EP), ' c_CaM4 = ', num2str(c_CaM4),...
             ' c_M = ', num2str(c_M), ' c_MP = ', num2str(c_MP), '\n']);
    error('Negative!');
end

% c_EP = c_E * paras.c_P / paras.K_EP;
% c_X = [c_CaM4, c_E, c_M, c_MP];

% Compute all the probabilities
if abs(c_Ms) < 1e-14
    p_M = 0; p_M_ATP = 0; p_M_CaM = 0; p_MA_CaM = 0; p_M_EP = 0;
else
    p_M = c_M / c_Ms;
    p_M_ATP = paras.c_ATP/paras.K_M_ATP * p_M;
    p_M_CaM = c_CaM4/paras.K_M_CaM * p_M;
    p_MA_CaM = c_CaM4/paras.K_MA_MACaM * p_M_ATP;
    p_M_EP = c_EP/paras.K_M_EP * p_M;
end

if abs(c_MPs) < 1e-14
    p_MP = 0; p_MP_ADP = 0; p_MP_CaM = 0; p_MPA_CaM = 0; p_MP_E = 0;
else
    p_MP = c_MP / c_MPs;
    p_MP_ADP = paras.c_ADP/paras.K_MP_ADP * p_MP;
    p_MP_CaM = c_CaM4/paras.K_MP_CaM * p_MP;
    p_MPA_CaM = c_CaM4/paras.K_MPA_MPACaM * p_MP_ADP;
    p_MP_E = c_E/paras.K_MP_E * p_MP;
end

end