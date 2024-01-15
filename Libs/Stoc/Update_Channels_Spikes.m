function [n_channels_new] = Update_Channels_Spikes(n_channels)
% AMPAR, NMDAR activation during pre-synaptic spikes
% INPUTS:
%        n_channels -- 1 a 35*1 array (number of channels different states)
% OUTPUTS:
%        n_channels_new -- 1 a 35*1 array

n_channels_new = n_channels;

n_AMPAR = n_channels(1:4); 
n_NMDAR = n_channels(5:8); % 1
[n_AMPAR_new] = Release_AMPAR(n_AMPAR);
[n_NMDAR_new] = Release_NMDAR(n_NMDAR);
n_channels_new(1:4) = n_AMPAR_new;
n_channels_new(5:8) = n_NMDAR_new;

end