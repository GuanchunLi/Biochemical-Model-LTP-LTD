function [n_channels] = Con2State_VC_II(y_Ca, paras)
%CON2NUM Convert a list of variables to numbers of channels
n_channels = zeros(22, 1);

% AMPAR number
s_AMPAR = y_Ca(2); % 1 Value of AMPAR gating variable s
n_AMPAR_x0 = mnrnd(paras.n_AMPAR, [1-s_AMPAR; s_AMPAR]);
n_AMPAR = [n_AMPAR_x0(1); 0; n_AMPAR_x0(2); 0];
n_channels(1:4) = n_AMPAR;

% NMDAR number
s_NMDAR = y_Ca(3); % 1 Value of NMDAR gating variable s
n_NMDAR_x0 = mnrnd(paras.n_NMDAR, [1-s_NMDAR; s_NMDAR]);
n_NMDAR = [n_NMDAR_x0(1); 0; n_NMDAR_x0(2); 0];
n_channels(5:8) = n_NMDAR;

% CaV12 number
y_CaV12 = y_Ca(4:9);  % Variables of CaV12 model
p_CaV12 = [1-sum(y_CaV12); y_CaV12(3); y_CaV12(2); y_CaV12(5); y_CaV12(4); y_CaV12(1); y_CaV12(6)];
n_CaV12 = mnrnd(paras.n_CaV12, p_CaV12);
n_channels(9:15) = n_CaV12;

% CaV13 number
y_CaV13 = y_Ca(10:15);  % Variables of CaV13 model
p_CaV13 = [1-sum(y_CaV13); y_CaV13(3); y_CaV13(2); y_CaV13(5); y_CaV13(4); y_CaV13(1); y_CaV13(6)];
n_CaV13 = mnrnd(paras.n_CaV13, p_CaV13);
n_channels(16:22) = n_CaV13;

% % Na channel number
% s_Na = y_Ca(18:19); % 1 Value of Sodium Current m, h
% m = s_Na(1); h = s_Na(2);
% m_dist = [(1-m)^3; 3*m*(1-m)^2; 3*m^2*(1-m); m^3];
% p_Na = [(1-h)*m_dist; h*m_dist];
% p_Na = p_Na / sum(p_Na);
% n_Na = mnrnd(paras.n_Na, p_Na);
% n_channels(23:30) = n_Na;

% % K channel number
% n = y_Ca(20); % 1 Value of Potassium Current n
% p_K = [(1-n)^4; 4*n*(1-n)^3; 6*(n^2)*(1-n)^2; 4*(1-n)*n^3; n^4];
% p_K = p_K / sum(p_K);
% n_K = mnrnd(paras.n_K, p_K);
% n_channels(31:35) = n_K;

end