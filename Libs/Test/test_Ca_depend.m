c_Ca_lst = 1e-2:1e-3:0.2;

K1 = 0.1; K2 = 0.5; k = 1; K_CaN = 0.1;
n = 4;
f_Ca = @(c_Ca) 4e2 * K_CaN ./ (1 + k * (1 + K1 ./ c_Ca).^n);
% f_Ca = @(c_Ca) 3e-10 * exp(150 * c_Ca);
semilogy(c_Ca_lst, f_Ca(c_Ca_lst));
hold on;

c_Ca_lst = [0.1, 0.11, 0.12, 0.14, 0.15, 0.18];
f_Ca_lst = [1e-4, 1e-2, 3e-2, 0.2, 0.8, 11];
semilogy(c_Ca_lst, f_Ca_lst, 'r*')


