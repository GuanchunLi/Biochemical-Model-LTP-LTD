clear;
addpath([pwd, '/../Parameters']);

paras.v = 1; 
paras = Param_CaM_Linse(paras);

% c_Ca_lst = sort([1e-4:1e-4:10, 10:1:1e3]);
c_Ca_lst = 5e-2:1e-4:1;
K_Ca_lst = 0*c_Ca_lst;
K_Ca_ori_lst = 0*c_Ca_lst;
for i = 1:length(c_Ca_lst)
    K_Ca_lst(i) = equal_constant_Ca(c_Ca_lst(i), paras);
    K_Ca_ori_lst(i) = equal_constant_Ca_test(c_Ca_lst(i), paras);
end
plot(c_Ca_lst, 1 - K_Ca_lst./(K_Ca_lst + 10^3));
hold on;
plot(c_Ca_lst, 1 - K_Ca_ori_lst./(K_Ca_ori_lst + 10^3));
plot(c_Ca_lst, (c_Ca_lst.^3)./(c_Ca_lst.^3 + 0.4^3));