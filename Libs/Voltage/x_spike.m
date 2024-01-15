function [x] = x_spike(t, tk_lst, tau, alpha)
tk_lst_eff = tk_lst(tk_lst < t);
x = alpha * sum(exp(-(t-tk_lst_eff) / tau));
end