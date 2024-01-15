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

%% Generate Calcium source

% parpool(32);

t_induce = [5, 40, 80, 115];
duration = [2, 10, 5, 3];
stim_mag = [0.3, 0.08, 0.1, 0.35];

[t_Ca_lst, Ca_lst, t_lst, y_lst, c_MP_lst, c_E_lst] = run_simulation_k(stim_mag, t_induce, duration, paras);
save('Result/Figures_For_Paper/Fig1/Model_dynamic_full.mat', 't_induce', 'duration', 'stim_mag', ...
 't_Ca_lst', 'Ca_lst', 't_lst', 'y_lst', 'c_MP_lst', 'c_E_lst');
fig = figure();
plot(t_lst, 100*c_MP_lst/(paras.n*paras.c0));
hold on;
plot(t_lst, 100*c_E_lst/paras.c_E_tot)
xlabel('Time (s) ', 'FontSize', 18);
ylabel('Activity (%)', 'FontSize', 18);
legend('CaMKII', 'PP');
% saveas(fig, 'Result/Result_04_28/Model/Model_dynamic_'+name_lst(k)+'.png');;


%% Aux functions

function [t_Ca_lst, Ca_lst, t_lst, y_lst, c_MP_lst, c_E_lst] = run_simulation_k(stim_mag, t_induce, duration, paras)
% t_induce = 25;
c_Ca_base = 0.1;
cal_source = @(t) c_Ca_base + stim_mag(1) * (t > t_induce(1)) .* (t < (t_induce(1)+duration(1))) ...
                            + stim_mag(2) * (t > t_induce(2)) .* (t < (t_induce(2)+duration(2))) ...
                            + stim_mag(3) * (t > t_induce(3)) .* (t < (t_induce(3)+duration(3))) ...
                            + stim_mag(4) * (t > t_induce(4)) .* (t < (t_induce(4)+duration(4)));
%% Initial state of CaMKII system
y0 = paras.y0_CaMKII;
% y0 = [0.155892302414455	90	99.8441076975867	4.97555338531664	12.8345512016881	5.06744168043888	0.0971183173176102	0.000418450032903543	0.000310254675171539	0.000155107349251540	1.76290921526365e-06	1.33675316639721e-06	1.33678274081244e-06	3.30423026229630e-07	7.42317787811958e-09	5.63239976647652e-09	2.88011338427329e-09	3.12634203774967e-11	2.82607553045392e-14	0.248354700099941]';
paras.c0 = sum(y0);

%% Start the simulation
y = y0; T = 160;
t_cri_lst = [t_induce(1), t_induce(1)+duration(1), t_induce(2), t_induce(2)+duration(2), ...
             t_induce(3), t_induce(3)+duration(3), t_induce(4), t_induce(4)+duration(4)];

[t_lst, y_lst] = Simu_ode_scheme_multistep(0, T, t_cri_lst, 1, 1e-4, y, paras, cal_source);

c_MP_lst = y_lst(:, 6:19) * paras.MP_lst';
c_E_lst= y_lst(:, 1);

t_Ca_lst = t_lst;
Ca_lst = cal_source(t_lst);

% y_ss = y_lst(end, :);
% CaMKII_SS = y_ss(6:19) * paras.MP_lst' / (paras.n * paras.c0) * 100;
% E_SS = y_ss(1); 
% AMPA_SS = y_ss(20) / y0(20);

% fprintf('%d finished! \n', k)
end
