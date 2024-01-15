function [state] = Con2State_CaMKII_family(paras, con)
%CON2NUM Convert a state of concentrations to a state of CaMKII units
state = zeros(paras.n+1, paras.n_CaMKII+1);

% PP number
con_E0 = con(1); con_EP0 = con(3);
num_PP = mnrnd(paras.n_PP, [con_E0, con_EP0]/(con_E0 + con_EP0));
num_E0 = num_PP(1); num_EP0 = num_PP(2);

% CaMKII number
con_CaMKII = con(5:19);
p_CaMKII = con_CaMKII / sum(con_CaMKII);
% sum(con0) * paras.coef_V must equal to paras.n_CaMKII!!
num_CaMKII = mnrnd(paras.n_CaMKII, p_CaMKII);

% Combine the states (see notes)

% Adding PP
state(1, end) = num_E0; state(2, end) = num_EP0;

% Adding CaMKII
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
k = 1;
for i = 1:length(num_CaMKII)
    if num_CaMKII(i) ~= 0
        state(:, k:(k+num_CaMKII(i)-1)) = state(:, k:(k+num_CaMKII(i)-1)) + S(:, i);
        k = k+num_CaMKII(i);
    end
end

%% AUX
% num0 = paras.coef_V * con0;
% num0_int = floor(num0);
% % Suppose sum(num0) = n_C0;
% k = paras.n_CaMKII - sum(num0_int);
% if k ~= 0
%     p = (num0 - num0_int)'; p = p / sum(p);
%     num0_k = mnrnd(k, p);
%     num0_final = num0_int + num0_k';
% else
%     num0_final = num0_int;
% end
end