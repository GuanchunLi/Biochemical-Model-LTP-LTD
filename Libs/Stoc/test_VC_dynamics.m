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
paras = Param_Channel_Numbers(paras);
paras.display = 0;

%%

% I_stim_func = @(t) 0.1;
tk_pre_lst = [1-10/1e3, 2+10/1e3];
tk_post_lst = [1, 2];
I_stim_func = @(t) 3000 * any((tk_post_lst > t+2/1e3) .* (tk_post_lst < t+3/1e3)) + 0 * (t>5);
NMDA_factor = @(t) 0.412;
CaV_factor = @(t) 2.02;

dt0 = 1e-4;

y_VC_0 = [-71.4479454884918; 0.103233668741414; 0.122470478619358; 0.00756800588549419; 0.980292454270282; 0.160041702948184];
n_channels_0 = Con2State_VC_II(paras.y0_Ca, paras);
% n_channels_0 = [20	0	0	0	9	0	0	0	0	0	0	0	0	27	3	0	0	0	0	0	45	3	1	0	0	0	48	4	0	0	5	7	0	0	0]';
% n_channels_0 = randi(10, 35, 1);

t_start = 0; t_end = 20;
% t = 1;
% [R_channels] = Rate_Channels(t, n_channels, y_VC(1), y_VC(3));
% [t_new, y_VC_new, n_channels_new, R_channels_new] = ...
%     Simu_VC_nl_one_step(t, y_VC, n_channels, R_channels, I_stim, dt, NMDA_factor, CaV_factor);

% [t_lst, y_VC_lst, n_channels_lst] = ...
%     Simu_VC_scheme_II(t_start, t_end, y_VC_0, n_channels_0, dt0, I_stim_func, NMDA_factor, CaV_factor);

[t_lst, y_VC_lst, n_channels_lst] = ...
    Simu_VC_scheme_Spikes_II(t_start, t_end, y_VC_0, n_channels_0, dt0, I_stim_func, NMDA_factor, CaV_factor, tk_pre_lst);