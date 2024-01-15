function K_Ca = equal_constant_Ca_test(c_Ca, paras)
    K_Ca = 1 ./ (c_Ca.^3) + 0.2 ./ (c_Ca.^4) + 1;
    % K_Ca = 100 ./ (c_Ca.^3) + 0.1 ./ (c_Ca.^3) + 0.2 ./ (c_Ca.^4) + 1;
end