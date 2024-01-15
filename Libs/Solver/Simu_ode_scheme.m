function [t_lst, y_lst, sol] = Simu_ode_scheme(t_start, t_end, max_step, y0, paras, cal_source)
% SIMU_ODE_SCHEME
% The scheme used to simulate the CaMKII reactions
options = odeset('RelTol',1e-7,'AbsTol',1e-9, 'MaxStep', max_step, 'NonNegative',1:length(y0), 'Stats', 'off');
% options = odeset('RelTol',1e-9,'AbsTol',1e-11, 'MaxStep', max_step, 'NonNegative',1:19, 'Stats', 'on');
% options = odeset('MaxStep', 1);
sol = ode15s(@(t, y) CaMKII_family_dynamic(t, y, cal_source(t), paras), [t_start, t_end], y0, options);
t_lst = sol.x;
y_lst = sol.y';
end