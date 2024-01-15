function [state0] = Con2State(paras, con0)
%CON2NUM Convert a state of concentrations to a state of CaMKII units 
num0 = paras.coef_V * con0;
num0_int = floor(num0);
% Suppose sum(num0) = n_C0;
k = paras.n_C0 - sum(num0_int);
if k ~= 0
    p = (num0 - num0_int)'; p = p / sum(p);
    num0_k = mnrnd(k, p);
    num0_final = num0_int + num0_k';
else
    num0_final = num0_int;
end
S = [0, 0, 0, 0, 0, 0, 0;
     1, 0, 0, 0, 0, 0, 0;
     1, 1, 0, 0, 0, 0, 0;
     1, 1, 1, 0, 0, 0, 0;
     1, 1, 0, 1, 0, 0, 0;
     1, 1, 0, 0, 1, 0, 0;
     1, 1, 1, 1, 0, 0, 0;
     1, 1, 1, 0, 1, 0, 0;
     1, 1, 1, 0, 0, 1, 0;
     1, 1, 0, 1, 0, 1, 0;
     1, 1, 1, 1, 1, 0, 0;
     1, 1, 1, 1, 0, 1, 0;
     1, 1, 1, 0, 1, 1, 0;
     1, 1, 1, 1, 1, 1, 0;
     1, 1, 1, 1, 1, 1, 1]';
state0 = zeros(paras.n+1, paras.n_C0);
k = 1;
for i = 1:length(num0_final)
    if num0_final(i) ~= 0
        state0(:, k:(k+num0_final(i)-1)) = state0(:, k:(k+num0_final(i)-1)) + S(:, i);
        k = k+num0_final(i);
    end
end

