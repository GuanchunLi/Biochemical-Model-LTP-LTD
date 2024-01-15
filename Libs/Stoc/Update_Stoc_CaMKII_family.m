function [t_new, y_new, Rate_new] = Update_Stoc_CaMKII_family(t, y, Rate_new_temp, T_mat, dt, Rate_func, paras)
% The reaction scheme for the CaMKII subunits with neighboring
% phosphorylation (Stochastic version)
% INPUTS:
%        t -- s Current time 
%        y -- a (1+6) * (n_CaMKII+1) matrix
%        c_Ca -- muM Concentration of Calcium at current time     
%        paras -- Group of parameters determing the fast equalibrium
% OUTPUTS:
%        R_mat_all -- a (1+5*6) * (1+n_CaMKII) matrix
% Dependency: 
%        compute_Prob_special.m

n = paras.n;
n_CaMKII = paras.n_CaMKII;
[tk, ik] = min(abs(T_mat), [], 1);
[tj, jk] = min(tk, [], 2);
idx_update = [ik(jk), jk];
dt_update = T_mat(idx_update(1), idx_update(2));

if dt_update > dt
    t_new = t + dt;
    y_new = y;
    Rate_new = Rate_new_temp;
elseif dt_update > 0
    % Update the reaction time
    t_new = t + dt_update;
    % Update y here to be y_new
    y_new = y;
    if idx_update(2) < (n_CaMKII+1) % CaMKII reactions
        if idx_update(1) == 1
            % Opening and closing of CaMKII
            y_new(idx_update(1), idx_update(2)) = 1 - y_new(idx_update(1), idx_update(2));
        else
            % Phosphorylation & Dephosphorylation of CaMKII
            k = mod(idx_update(1)-1, n);
            if k == 0
                k = k + paras.n;
            end
            y_new(k+1, idx_update(2)) = 1 - y_new(k+1, idx_update(2));
            if idx_update(1) >= (5*n+2)
            % Dephosphorylatio of CaMKII by PP
                if y(k+1, idx_update(2)) == 1
                    y_new(1, end) = y_new(1, end) - 1;
                    y_new(2, end) = y_new(2, end) + 1;
                else
                    y_new(1, end) = y_new(1, end) + 1;
                    y_new(2, end) = y_new(2, end) - 1;
                end
            end
        end
    else % PP reactions
        if idx_update(1) <= 2
            % Auto/Calcium dephosphorylation of PP
            y_new(1, end) = y_new(1, end) + 1;
            y_new(2, end) = y_new(2, end) - 1;
        else
            error("PP Reaction error!");
        end
    end
    % Update the reaction rate
    Rate_new = Rate_func(t_new, y_new);
else
    error("dt < 0!");
end

end