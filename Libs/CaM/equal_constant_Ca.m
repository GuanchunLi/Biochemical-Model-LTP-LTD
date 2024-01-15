function K_Ca = equal_constant_Ca(c_Ca, paras)
    K_1 = paras.K_1; K_2 = paras.K_2;
    K_3 = paras.K_3; K_4 = paras.K_4;
    c_4 = 1; c_3 = K_4 / c_Ca; c_2 = c_3 * K_3 / c_Ca;
    c_1 = c_2 * K_2 / c_Ca; c_0 = c_1 * K_1 / c_Ca;
    K_Ca = c_0 + c_1 + c_2 + c_3 + c_4;
end