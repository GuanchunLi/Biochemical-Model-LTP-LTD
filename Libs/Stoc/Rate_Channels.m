function [R_channels] = Rate_Channels(t, n_channels, V, cp)
% The reaction scheme for the CaMKII subunits with neighboring
% phosphorylation (Stochastic version)
% INPUTS:
%        t -- s Current time 
%        n_channels -- 1 a 35*1 array
%        V -- mC membrane potential
%        cp -- muM Concentration of Calcium (of subspace) at current time     
% OUTPUTS:
%        R_channels -- a 84*1 array

n_AMPAR = n_channels(1:4); % 1 
n_NMDAR = n_channels(5:8); % 1
n_CaV12 = n_channels(9:15); % 1
n_CaV13 = n_channels(16:22); % 1
n_Na = n_channels(23:30); % 1
n_K = n_channels(31:35); % 1

R_AMPAR = Rate_AMPAR(n_AMPAR);
R_NMDAR = Rate_NMDAR(n_NMDAR);
R_CaV12 = Rate_CaV12_channel(n_CaV12, V, cp);
R_CaV13 = Rate_CaV13_channel(n_CaV13, V, cp);
R_Na = Rate_Na_channel(n_Na, V);
R_K = Rate_K_channel(n_K, V);

R_channels = zeros(78, 1);
R_channels(1:5) = R_AMPAR;
R_channels(6:10) = R_NMDAR;
R_channels(11:30) = R_CaV12;
R_channels(31:50) = R_CaV13;
R_channels(51:70) = R_Na;
R_channels(71:78) = R_K;

R_channels(R_channels < 0) = 0;

end