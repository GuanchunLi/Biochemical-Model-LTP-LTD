a4 = 0.2; a3 = 1; a2 = 10; a1 = 10; 

paras.K_1 = a4 / a3;
paras.K_2 = a3 / a2;
paras.K_3 = a2 / a1;
paras.K_4 = a1;

F= @(X) [X(2) + X(4) - a1;
         X(1)*X(2) + X(2)*X(4) + X(3)*X(4) - a2;
         (X(1)+X(3)) * X(2)*X(4) - a3/a1;
         X(1)*X(2)*X(3)*X(4) - a4/a2];

x = lsqnonlin(F, [0.1, 0.2, 2, 2], [0, 0, 0, 0], [inf, inf, inf, inf]);
K_1_C = x(1)
K_2_C = x(2)
K_1_N = x(3)
K_2_N = x(4)