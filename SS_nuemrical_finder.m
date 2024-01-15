
%% Initialization & Parameter setting
clear;

addpath(genpath([pwd, '/Libs/']));
addpath(genpath([pwd, '/Parameters']));

%%
% Set up all parameters for the model
paras.n = 6;    % Number of monomers in each CaMKII (consider only one ring)
paras.v = 1;    % Reaction speed scale (reaction_rate = base_rate * v)
paras.v_tot = 2;
paras = Param_CaM(paras);
paras = Param_Concent_Retuning(paras);
paras = Param_FastEqual(paras);
paras = Param_ReactRate_Retuning(paras);
paras = Param_PhosTrans(paras);
paras = Param_PP1_mod3_Retuning(paras);
paras.display = 0;
paras.c0 = 18;
paras = Param_RateModel(paras);

% Define the time course of calcium source
c_Ca = 0.1;
c_E_tot = 100; %100;
% c_E0 = 0.1; c_EP0 = c_E_tot - c_E0;

% c_E0_lst = sort([1e-3:2e-3:1e-1, 1e-1:2e-1:20, 20:2:100]);

c_E0_lst = unique(sort([1e-2:1e-2:2, 1:1e-1:10, 10:1:100]));
% c_E0_lst = [0:1e-2:2];
% c_E0_lst = [2];

%% Another half of the phase plane
% c_E0_lst = sort([7.9e-2:1e-5:8e-2]);
c_E0_lst = unique(sort([7.9e-2:1e-5:8e-2, 1e-2:1e-2:2, 1:1e-1:10, 10:1:99, 99:1e-1:100]));
c_MP_lst = 0 * c_E0_lst;
for i = 1:length(c_E0_lst)
    c_E0 = c_E0_lst(i);
    c_EP0 = c_E_tot - c_E0;
    y_E_SS = [c_E0; 100; c_EP0; 5];
    [c_MP] = E_SS_solver_new_mod3_Retuning(y_E_SS, paras, c_Ca);
    c_MP_lst(i) = c_MP;
end
% figure();
semilogx(c_MP_lst, c_E0_lst, 'k-');
hold on;
% xlim([1e-9, 1e-1]);
% ylim([0.945, 0.95]);


%% TEST CODE
% y0 = zeros(15, 1);
% y0(1, :) = 12*0.93;
% y0(2, :) = 12*0.01;
% y0(3, :) = 12*0.06;
% c_E0 = 0.949244; c_EP0 = c_E_tot - c_E0;
% y_CaMKII_SS = CaMIKK_SS_finder(y0, c_Ca, c_E0, c_E_tot - c_E0, paras);
% [K_MP_E_1, K_MP_E_2, p_MP] ...
%     = ReadOut_CaMKII_System(100, y_CaMKII_SS, c_Ca, c_E0, c_EP0, paras.c_CaM, paras);
% c_MPs = y_CaMKII_SS(2:end) * paras.MP_lst';
% c_MP = p_MP * c_MPs
% 
% y_E_SS = [c_E0; 90; c_EP0; 5];
% [c_MP_2] = E_SS_solver(y_E_SS, paras)

%% Half of the phase plane
% Set initial value
y0_low = zeros(15, 1);
y0_low(1, :) = paras.c0*0.93; y0_low(3, :) = paras.c0*0.01; y0_low(4, :) = paras.c0*0.06;

y0_high = zeros(15, 1);
y0_high(15, :) = paras.c0*0.93; y0_high(14, :) = paras.c0*0.01; y0_high(13, :) = paras.c0*0.06;

y0_low = y0_low; y0_high = y0_high;
% y0_low = [11.6751036611599,0.0497836063147355,0.100115745785740,0.0454417265424090,0.0338800733909399,0.0167945677703288,0.0202931806016947,0.0153690179695917,0.0153848308435875,0.00385824164946707,0.00906180295434331,0.00692584681319082,0.00351618326429362,0.00407773082075718,0.000393790552688682]';
% y0_high = [7.97478440599491,0.00160596091048060,0.0147972497604925,0.0463224184765831,0.0258647958674598,0.0120384286582611,0.142019728962261,0.0820402994933283,0.0824814611126153,0.0167705213663273,0.441417522267031,0.274291225964399,0.139673993588341,1.46618593788761,1.28079239352621]';

c_MPs_lst_low = 0*c_E0_lst;
c_Ms_lst_low = 0*c_E0_lst;
c_MC_lst_low = 0*c_E0_lst;
c_MP_lst_low = 0*c_E0_lst;

c_MPs_lst_high = 0*c_E0_lst;
c_Ms_lst_high = 0*c_E0_lst;
c_MC_lst_high = 0*c_E0_lst;
c_MP_lst_high = 0*c_E0_lst;
% 

y_CaMKII_SS_low = {};
y_CaMKII_SS_high = {};

parfor i = 1:length(c_E0_lst)
    c_E0 = c_E0_lst(i);
    c_EP0 = c_E_tot - c_E0;
    y_CaMKII_SS = CaMIKK_SS_finder(y0_low, c_Ca, c_E0, c_EP0, paras);
    % fprintf('C=%f, O0=%f \n', y_CaMKII_SS(1), y_CaMKII_SS(2));
    [K_MP_E_1, K_MP_E_2, p_MP] ...
    = ReadOut_CaMKII_System(100, y_CaMKII_SS, c_Ca, c_E0, c_EP0, paras.c_CaM, paras);
    c_MPs = y_CaMKII_SS(2:end) * paras.MP_lst';
    c_Ms = y_CaMKII_SS(2:end) * paras.M_lst';
    c_MC = y_CaMKII_SS(1);
    c_MPs_lst_low(i) = c_MPs; c_Ms_lst_low(i) = c_Ms; c_MC_lst_low(i) = c_MC;
    c_MP_lst_low(i) = p_MP * c_MPs;
    y_CaMKII_SS_low{i} = y_CaMKII_SS;

    y_CaMKII_SS = CaMIKK_SS_finder(y0_high, c_Ca, c_E0, c_EP0, paras);
    % fprintf('C=%f, O0=%f \n', y_CaMKII_SS(1), y_CaMKII_SS(2));
    [K_MP_E_1, K_MP_E_2, p_MP] ...
    = ReadOut_CaMKII_System(100, y_CaMKII_SS, c_Ca, c_E0, c_EP0, paras.c_CaM, paras);
    c_MPs = y_CaMKII_SS(2:end) * paras.MP_lst';
    c_Ms = y_CaMKII_SS(2:end) * paras.M_lst';
    c_MC = y_CaMKII_SS(1);
    c_MPs_lst_high(i) = c_MPs; c_Ms_lst_high(i) = c_Ms; c_MC_lst_high(i) = c_MC;
    c_MP_lst_high(i) = p_MP * c_MPs;
    y_CaMKII_SS_high{i} = y_CaMKII_SS;
end

fprintf('Finished! \n');
% figure();
c_MP_lst_tot = [c_MP_lst_low, c_MP_lst_high];
c_E0_lst_tot = [c_E0_lst, c_E0_lst];
[c_MP_lst_tot, sort_ind] = sort(c_MP_lst_tot);
c_E0_lst_tot = c_E0_lst_tot(sort_ind);
% plot(c_MP_lst_low, c_E0_lst, 'bo');
% plot(c_MP_lst_high, c_E0_lst, 'bo');
semilogx(c_MP_lst_tot, c_E0_lst_tot, 'b-');
hold on;
xlabel('c_{MP}');
ylabel('[E_0]');
% xlim([0, 1e-4]);
% plot(c_MPs_lst_low, c_E0_lst, 'bo');
% plot(p_MPs_lst_high, c_E0_lst, 'bo');
% xlim([5.9e-3, 5.95e-3]);

% figure();
% c_MPs_lst_tot = [c_MPs_lst_low, c_MPs_lst_high];
% c_E0_lst_tot = [c_E0_lst, c_E0_lst];
% [c_MPs_lst_tot, sort_ind] = sort(c_MPs_lst_tot);
% c_E0_lst_tot = c_E0_lst_tot(sort_ind);
% plot(c_E0_lst_tot, c_MPs_lst_tot, 'k.');

%% Aux functions

function [y_CaMKII_SS, t_lst, y_lst] = CaMIKK_SS_finder(y_CaMKII_0, c_Ca, c_E0, c_EP0, paras)
    T = 3000;
    Reaction_scheme = @(t, y_CaMKII) ...
        Reaction_ode_wi_transpp_RE(t, y_CaMKII, c_Ca, c_E0, c_EP0, paras.c_CaM, paras);
    options = odeset('RelTol',1e-9,'AbsTol',1e-8, 'NonNegative',1:15, 'Stats', 'off');
    sol = ode15s(Reaction_scheme, [0, T], y_CaMKII_0, options);
    t_lst = sol.x;
    y_lst = sol.y';
    y_CaMKII_SS = y_lst(end, :);
end

function [dy_PP1] = PP1_dynamic_fixing(t, y_PP1,  c_Ca, c_Ms, c_MPs, paras)
    c_E0 = y_PP1(1);
    c_EP0 = y_PP1(3);
    [p_M, p_M_ATP, p_M_CaM, p_MA_CaM, p_M_EP, ...
     p_MP, p_MP_ADP, p_MP_CaM, p_MPA_CaM, p_MP_E, c_CaM4, c_E, c_EP] ...
     = compute_Prob(c_Ms, c_MPs, c_Ca, c_E0, c_EP0, paras.c_CaM, paras);
    p_E = c_E  / c_E0;
    p_EP = c_EP / c_EP0;
    gamma_c_MPs = (c_MPs + paras.c_MPs_base) / c_MPs;
    
    dE_CaM = p_M_EP * paras.beta_E * c_Ms - p_MP_E * paras.alpha_E * c_MPs * gamma_c_MPs;
    % Compute the avaliable active CaM
    K_Ca = equal_constant_Ca(c_Ca, paras.K_1, paras.K_2, paras.K_3, paras.K_4);
    c_CaM4_tot = paras.c_CaM - (K_Ca - 1) * c_CaM4;
    dy_PP1 = PP1_dynamic_new(t, y_PP1, p_E, p_EP, dE_CaM, c_CaM4_tot, paras);
end

function [y_E_SS, t_lst, y_lst] = E_SS_finder(y_E_0, c_Ca, c_Ms, c_MPs, paras)
    T = 50;
    Reaction_scheme = @(t, y_PP1) ...
        PP1_dynamic_fixing(t, y_PP1,  c_Ca, c_Ms, c_MPs, paras);
    options = odeset('RelTol',1e-3,'AbsTol',1e-7, 'NonNegative',1:4, 'Stats', 'off');
    sol = ode15s(Reaction_scheme, [0, T], y_E_0, options);
    t_lst = sol.x;
    y_lst = sol.y';
    y_E_SS = y_lst(end, :);
end


function [c_MP] = E_SS_solver(y_E_SS, paras, c_Ca)
    c_E0 = y_E_SS(1);
    c_EP0 = y_E_SS(3);
    alpha_P_auto_Ca = 0.1 * 10; % alpha_P_auto_func(c_Ca, paras);
    alpha_P_Ca = 0.02 * 0.1; % alpha_P_func(c_Ca, paras);
%     c_E = (c_E0 - (paras.d_E_base)/paras.alpha_E) / (1 + ...
%         paras.alpha_P/paras.alpha_E * c_EP0 / (paras.K_m50_EP + c_EP0));
    c_E = (paras.alpha_E*c_E0 - paras.d_E_base - alpha_P_Ca*c_EP0 - 10*c_E0) / ...
        (paras.alpha_E + alpha_P_auto_Ca*c_EP0/(paras.K_m50_EP + c_EP0) - paras.r_P_base);
    c_MP_E = c_E0 - c_E;
    c_MP = paras.K_MP_E * c_MP_E / c_E;
end

function [c_MP] = E_SS_solver_new_mod1(y_E_SS, paras, c_Ca)
% % paras.alpha_P_auto = 10; % 1000 1/s PP1 auto-dephosphorylation rate
% % paras.alpha_P = 0; % 1/s PP1 dephosphorylation rate through Ca pathway
% % 
% % paras.alpha_Cdk5 = 220; % 1/s Cdk5-PP1 phosphorylation rate
% % paras.K_m50_EP = 1e-4; % muM PP1 auto-dephosphorylation Michalis Menton Coef
% % paras.K_m50_E = 20; % muM Cdk5-PP1 phosphorylation Michalis Menton Coef

    c_E0 = y_E_SS(1);
    c_EP0 = y_E_SS(3);
    K_m50_EP = paras.K_m50_EP;
    K_m50_E = paras.K_m50_E;
    % alpha_P_Ca = alpha_P_func(c_Ca, paras);
    Cdk5_a_Ca = Cdk5_a_func(c_Ca, paras);
    A = paras.alpha_P_auto * c_EP0 / (K_m50_EP + c_EP0) + paras.alpha_E;
    B = paras.alpha_E * c_E0 - paras.d_E_base;
    C = paras.alpha_Cdk5 * Cdk5_a_Ca;
    % c_E = (B + C * (c_E0/(c_E0+K_m50_E))) / A;
    D = A * K_m50_E - B - C;
    c_E = (-D + sqrt(D^2 + 4*A*B*K_m50_E)) / (2*A);
    if c_E < 0
        c_E = c_E0 * (1e-10);
    elseif c_E >= c_E0
        c_E = c_E0 * (1-1e-10);
    end
    c_MP_E = c_E0 - c_E;
    c_MP = paras.K_MP_E * c_MP_E / c_E;
end

function [c_MP] = E_SS_solver_new_mod2(y_E_SS, paras, c_Ca)
% %     paras.alpha_P_auto = 10; % 1000 1/s PP1 auto-dephosphorylation rate
% %     paras.alpha_P = 0; % 1/s PP1 dephosphorylation rate through Ca pathway
% % 
% %     paras.alpha_Cdk5 = 1 * (10 + 5e-3) * 1e4; % 1/s Cdk5-PP1 phosphorylation rate
% %     paras.K_m50_EP = 1e-4; % muM PP1 auto-dephosphorylation Michalis Menton Coef
% %     paras.K_m50_E = 1e4;
    c_E0 = y_E_SS(1);
    c_EP0 = y_E_SS(3);
    K_m50_EP = paras.K_m50_EP;
    K_m50_E = paras.K_m50_E;
    % alpha_P_Ca = alpha_P_func(c_Ca, paras);
    Cdk5_a_Ca = Cdk5_a_func(c_Ca, paras);
    A = paras.alpha_P_auto * c_EP0 / (K_m50_EP + c_EP0) + paras.alpha_E;
    B = paras.alpha_E * c_E0 - paras.d_E_base;
    C = paras.alpha_Cdk5 * Cdk5_a_Ca;
    c_E = (B + C * (c_E0/(c_E0+K_m50_E))) / A;
    if c_E < 0
        c_E = c_E0 * (1e-10);
    elseif c_E >= c_E0
        c_E = c_E0 * (1-1e-10);
    end
    c_MP_E = c_E0 - c_E;
    c_MP = paras.K_MP_E * c_MP_E / c_E;
end

function [c_MP] = E_SS_solver_new_mod3(y_E_SS, paras, c_Ca)
    c_E0 = y_E_SS(1);
    c_EP0 = y_E_SS(3);
    K_m50_EP = paras.K_m50_EP;
    alpha_P_Ca = alpha_P_func(c_Ca, paras);
    alpha_E = paras.alpha_E * (1 + paras.c_MP_base/paras.K_MP_E);
    c_MP_E = (paras.alpha_P_auto * c_E0 * c_EP0 / (K_m50_EP + c_EP0) ...
              + paras.d_E_base + alpha_P_Ca * c_EP0) ...
           / (alpha_E + paras.alpha_P_auto * c_EP0 / (K_m50_EP + c_EP0));
    c_E = c_E0 - c_MP_E;
    c_MP = paras.K_MP_E * c_MP_E / c_E;
end

function [c_MP] = E_SS_solver_new_mod3_Retuning(y_E_SS, paras, c_Ca)
    c_E0 = y_E_SS(1);
    c_EP0 = y_E_SS(3);
    K_m50_EP = paras.K_m50_EP;
    alpha_P_Ca = alpha_P_func_Retuning(c_Ca, paras);
    alpha_E = paras.alpha_E * (1 + paras.c_MP_base/paras.K_MP_E);
    c_MP_E = (paras.alpha_P_auto * c_E0 * c_EP0 / (K_m50_EP + c_EP0) ...
              + paras.d_E_base + alpha_P_Ca * c_EP0) ...
           / (alpha_E + paras.alpha_P_auto * c_EP0 / (K_m50_EP + c_EP0));
    c_E = c_E0 - c_MP_E;
    c_MP = paras.K_MP_E * c_MP_E / c_E;
end

% function [c_MP] = E_SS_solver_new_mod3(y_E_SS, paras, c_Ca)
%     c_E0 = y_E_SS(1);
%     c_EP0 = y_E_SS(3);
%     K_m50_EP = paras.K_m50_EP;
%     alpha_P_Ca = alpha_P_func(c_Ca, paras);
%     c_E_func = @(c_E) ...
%     (paras.alpha_P_auto * c_EP0 / (K_m50_EP + c_EP0) ...
%      - paras.alpha_E*paras.c_MP_base/paras.K_MP_E) * c_E ...
%       + alpha_P_Ca * c_EP0 - paras.alpha_E*(c_E0-c_E) + paras.d_E_base;
%     if c_E_func(c_E0) * c_E_func(0) < 0
%         c_E = fzero(c_E_func, [0, c_E0]);
%     elseif c_E_func(c_E0) > 0
%         c_E = c_E0 * (1-1e-10);
%     elseif c_E_func(c_E0) <= 0
%         c_E = c_E0 * (1e-10);
%     end
%     c_MP_E = c_E0 - c_E;
%     c_MP = paras.K_MP_E * c_MP_E / c_E;
% end


function [c_MP] = E_SS_solver_new_mod6(y_E_SS, paras, c_Ca)
    c_E0 = y_E_SS(1);
    c_EP0 = y_E_SS(3);
    K_m50_EP = paras.K_m50_EP;
    alpha_P_Ca = alpha_P_func(c_Ca, paras);
    alpha_E = paras.alpha_E * (1 + paras.c_MP_base/paras.K_MP_E);
    c_MP_E = (paras.alpha_P_auto * c_E0 * c_EP0 / (K_m50_EP + c_EP0) ...
              + paras.d_E_base + alpha_P_Ca * c_EP0 - paras.r_P_base*c_E0) ...
           / (alpha_E + paras.alpha_P_auto * c_EP0 / (K_m50_EP + c_EP0));
    c_E = c_E0 - c_MP_E;
    c_MP = paras.K_MP_E * c_MP_E / c_E;
end

function [c_MP] = E_SS_solver_new_mod7(y_E_SS, paras, c_Ca)
    c_E0 = y_E_SS(1);
    c_EP0 = y_E_SS(3);
    K_m50_EP = K_m50_EP_func(c_Ca, paras);
    alpha_P_Ca = alpha_P_func_mod7(c_Ca, paras);
    c_MP_E = (paras.alpha_P_auto * c_E0 * c_EP0 / (K_m50_EP + c_EP0) ...
              + paras.d_E_base + alpha_P_Ca * c_EP0) ...
           / (paras.alpha_E + paras.alpha_P_auto * c_EP0 / (K_m50_EP + c_EP0));
    c_E = c_E0 - c_MP_E;
    c_MP = paras.K_MP_E * c_MP_E / c_E;
end

function [fitresult, xData, yData] = fit_SS_curve(x_lst, y_lst)
[xData, yData] = prepareCurveData(x_lst, y_lst );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 1;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

end
