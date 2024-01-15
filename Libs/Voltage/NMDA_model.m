function [dy] = NMDA_model(y, Vs)

x_NMDA = 1;
s_NMDA = y;

ds_NMDA = NMDA(s_NMDA, x_NMDA);
% [I_NMDA, I_NMDA_Ca] = NMDA_current(s_NMDA, Vs);
% I_NMDA = I_NMDA * NMDA_factor;
% I_NMDA_Ca = I_NMDA_Ca * NMDA_factor;

dy = ds_NMDA;
end