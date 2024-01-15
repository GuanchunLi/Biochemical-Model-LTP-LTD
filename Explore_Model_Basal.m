%% Initialization
clear;

addpath(genpath([pwd, '/Libs/']));
addpath(genpath([pwd, '/Parameters']));

% Set up all parameters for the model
paras.n = 6;    % Number of monomers in each CaMKII (consider only one ring)
paras.v = 1;    % Reaction speed scale (reaction_rate = base_rate * v)
% paras = Param_CaM(paras);
paras.v_tot = 2;
paras = Param_CaM(paras);
paras = Param_Concent_Retuning(paras);
paras = Param_FastEqual(paras);
paras = Param_ReactRate_Retuning(paras);
paras = Param_PhosTrans(paras);
paras = Param_PP1_mod3_Retuning(paras);
paras = Param_RateModel(paras);
paras = Param_Readout(paras);
paras = Param_InitVal_Retuning(paras);
paras.c0 = 18;
paras.display = 0;

% paras.alpha_P = 1.2;
% paras.alpha_P_auto = 0.;
%% Generate Calcium source

% parpool(32);

duration_lst = [10];
mag_lst = [0.05];
% duration_lst = [20];
% mag_lst = [1];
name_lst = ["LTD_stable"];
% N_trail = length(mag_lst);
% CaMKII_lst = zeros(N_trail, 1);
% E_lst = zeros(N_trail, 1);
% AMPA_lst = zeros(N_trail, 1);

for k = 1:length(duration_lst)
    duration = duration_lst(k);
    stim_mag = mag_lst(k);
    [t_Ca_lst, Ca_lst, t_lst, y_lst, c_MP_lst, c_Ms_lst, c_E_lst] = run_simulation_k(stim_mag, duration, paras);
    save('Result/Result_04_28/Model/Model_dynamic_'+name_lst(k)+'.mat', 'duration', 'stim_mag', ...
     't_Ca_lst', 'Ca_lst', 't_lst', 'y_lst', 'c_MP_lst', 'c_E_lst');
    fig = figure();
    plot(t_lst, c_MP_lst);
    hold on;
    plot(t_lst, c_E_lst)
    
    xlabel('Time (s) ', 'FontSize', 18);
    ylabel('Activity (%)', 'FontSize', 18);
    legend('CaMKII', 'PP');
    % saveas(fig, 'Result/Result_04_28/Model/Model_dynamic_'+name_lst(k)+'.png');
    fprintf('%d finished! \n', k);

%     figure();
%     plot(t_lst, [c_Ms_lst, c_MP_lst]);
end

%% Aux functions

function [t_Ca_lst, Ca_lst, t_lst, y_lst, c_MP_lst, c_Ms_lst, c_E_lst] = run_simulation_k(stim_mag, duration, paras)
t_induce = 10;
c_Ca_base = 0.1;
cal_source = @(t) c_Ca_base + stim_mag * (t > t_induce) .* (t < (t_induce+duration));
%% Initial state of CaMKII system

% y0 = paras.y0_CaMKII;
% y0 = [0.155892302414455	90	99.8441076975867	4.97555338531664	12.8345512016881	5.06744168043888	0.0971183173176102	0.000418450032903543	0.000310254675171539	0.000155107349251540	1.76290921526365e-06	1.33675316639721e-06	1.33678274081244e-06	3.30423026229630e-07	7.42317787811958e-09	5.63239976647652e-09	2.88011338427329e-09	3.12634203774967e-11	2.82607553045392e-14	0.248354700099941]';
% LTP
% y0 = [0.0799238593236288	90	99.9200761406775	4.97555338531664	7.93541731472326	2.37337195956760	1.05883727920599	1.16948743063646	0.270171419178263	0.0949232981807396	1.09494021836432	0.297268900645660	0.321791571014307	0.0381981158770563	1.03571332621252	0.392587231848421	0.215621186130208	1.18704726260487	0.514622901257440	0.315567556519906]';
% LTD
y0 = [99.9097289718654	90	0.0902710281358983	4.97555338531664	12.4806740543511	4.89324576774183	0.601079098536288	0.00614113043851307	0.0123268088638881	0.00615283224255007	4.19501466233275e-05	0.000126080545203678	0.000125937436026351	8.44078744380928e-05	2.16227976942845e-07	8.63267438632255e-07	6.45730085555501e-07	9.04485985164368e-10	4.45479388171807e-15	0.0421048203308460]';
paras.c0 = sum(y0);

%% Start the simulation
y = y0; T = 40;
t_cri_lst = [t_induce, t_induce+duration];
[t_lst, y_lst] = Simu_ode_scheme_multistep(0, T, t_cri_lst, 1, 1e-4, y, paras, cal_source);

c_MP_lst = y_lst(:, 6:19) * paras.MP_lst';
c_Ms_lst = y_lst(:, 6:19) * paras.M_lst';
c_E_lst= y_lst(:, 1);

t_Ca_lst = t_lst;
Ca_lst = cal_source(t_lst);

% y_ss = y_lst(end, :);
% CaMKII_SS = y_ss(6:19) * paras.MP_lst' / (paras.n * paras.c0) * 100;
% E_SS = y_ss(1); 
% AMPA_SS = y_ss(20) / y0(20);

% fprintf('%d finished! \n', k)
end
