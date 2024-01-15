function [t_new, y_new, Rate_new] = Simu_stoch_nl_one_step(t, y, Rate, Rate_func, Update_func, dt)
% One step of stochastic simulation of nonlinear poisson process with
% time-varying rate
U = rand(size(Rate));
Rate_new_temp = Rate_func(t+dt, y);
Delta = Rate.^2 - 2 * (Rate_new_temp-Rate) .* log(U) / dt;
T = 0 * Rate;
index = Delta > 0;
T(~index) = 2*dt;
T(index) = - 2 * log(U(index)) ./ (sqrt(Delta(index)) + Rate(index));
[t_new, y_new, Rate_new] = Update_func(t, y, Rate_new_temp, T, dt, Rate_func);
end