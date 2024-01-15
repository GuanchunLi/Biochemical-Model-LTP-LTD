function [n_channels_new] = Update_Channels(n_channels, i_channels)
% The reaction scheme for the CaMKII subunits with neighboring
% phosphorylation (Stochastic version)
% INPUTS:
%        n_channels -- 1 a 35*1 array (number of channels different states)
%        i_channels -- 1 index of the next reaction (1<=i_channels<=78)
% OUTPUTS:
%        n_channels_new -- 1 a 35*1 array

n_channels_new = n_channels;

if i_channels <= 5 % AMPAR reaction
    n_AMPAR = n_channels(1:4); 
    i_AMPAR = i_channels;
    [n_AMPAR_new] = Update_AMPAR(n_AMPAR, i_AMPAR);
    n_channels_new(1:4) = n_AMPAR_new;

elseif i_channels <= 10 % NMDAR reaction
    n_NMDAR = n_channels(5:8); % 1
    i_NMDAR = i_channels - 5;
    [n_NMDAR_new] = Update_NMDAR(n_NMDAR, i_NMDAR);
    n_channels_new(5:8) = n_NMDAR_new;
   
elseif i_channels <= 30 % CaV12 reaction
    n_CaV12 = n_channels(9:15); % 1
    i_CaV12 = i_channels - 10;
    [n_CaV12_new] = Update_CaV12_channel(n_CaV12, i_CaV12);
    n_channels_new(9:15) = n_CaV12_new;

elseif i_channels <= 50 % CaV13 reaction
    n_CaV13 = n_channels(16:22); % 1
    i_CaV13 = i_channels - 30;
    [n_CaV13_new] = Update_CaV13_channel(n_CaV13, i_CaV13);
    n_channels_new(16:22) = n_CaV13_new;

elseif i_channels <= 70 % Na channel reaction
    n_Na = n_channels(23:30); % 1
    i_Na = i_channels - 50;
    [n_Na_new] = Update_Na_channel(n_Na, i_Na);
    n_channels_new(23:30) = n_Na_new;

elseif i_channels <= 78 % K channel reaction
    n_K = n_channels(31:35); % 1
    i_K = i_channels - 70;
    [n_K_new] = Update_K_channel(n_K, i_K);
    n_channels_new(31:35) = n_K_new;
else
    error('Not valid reaction index!');
end
end