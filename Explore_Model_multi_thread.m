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
paras.display = 0;

%% Generate Calcium source

parpool(32);

duration_lst = 0.1:0.01:5;
mag_lst = 0.01:0.0005:1;
CaMKII_mat = zeros(length(duration_lst), length(mag_lst));
E_mat = 0*CaMKII_mat;
AMPA_mat = 0*CaMKII_mat;
n_index = 1:length(mag_lst);

parfor k = 1:length(duration_lst)
    for n = n_index
        duration = duration_lst(k);
        stim_mag = mag_lst(n);
        [CaMKII_SS, E_SS, AMPA_SS] = run_simulation_k(stim_mag, duration, paras);
        CaMKII_mat(k, n) = CaMKII_SS;
        E_mat(k, n) = E_SS;
        AMPA_mat(k, n) = AMPA_SS;
    end
    fprintf('%d finished! \n', k);
end

fig1 = figure(1);
heatmap(AMPA_mat);
saveas(fig1, 'Pulse_Readout.png');

save('Pulse_Readout.mat', 'duration_lst', 'mag_lst', 'CaMKII_mat', 'E_mat', 'AMPA_mat');

%% Aux functions

function [CaMKII_SS, E_SS, AMPA_SS] = run_simulation_k(stim_mag, duration, paras)
t_induce = 20;
cal_source = @(t) 0.1 + stim_mag * (t > t_induce) * (t < (t_induce+duration));

%% Initial state of CaMKII system
y0 = paras.y0_CaMKII;

%% Start the simulation
y = y0; T = 200;
t_cri_lst = [t_induce, t_induce+duration];
[t_lst, y_lst] = Simu_ode_scheme_multistep(0, T, t_cri_lst, 1, 1e-4, y, paras, cal_source);

y_ss = y_lst(end, :);
CaMKII_SS = y_ss(6:19) * paras.MP_lst' / (paras.n * paras.c0) * 100;
E_SS = y_ss(1); 
AMPA_SS = y_ss(20) / y0(20);

% fprintf('%d finished! \n', k)
end
