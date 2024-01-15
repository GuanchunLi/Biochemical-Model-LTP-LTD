clear;
y0_Ca = [-71.4479454884918	4.90000000000000e-323	2.17000000000000e-321	0.941129111726702	0.000124442334267052	2.53974631380274e-06	3.97962603451610e-06	0.000229440553805646	0.0584980417794499	0.527243225961749	0.00296438261162864	6.05001493656064e-05	3.15617570090125e-05	0.00546558044018662	0.463938310818898	0.103233668741414	0.122470478619358	0.00756800588549419	0.980292454270282	0.160041702948184]';

Vs = y0_Ca(1); % mV Membrane voltage of spine head
% s_AMPA = y_cir(2); % 1 Value of AMPA gating variable s
s_NMDA = y0_Ca(3); % 1 Value of NMDA gating variable s
y0_CaV12 = y0_Ca(4:9);  % Variables of CaV12 model
y0_CaV13 = y0_Ca(10:15); % Variables of CaV13 model
c_Ca = y0_Ca(16); % muM Intracelluar Calcium concentration
c_Ca_SS = y0_Ca(17); % muM Subspace Calcium concentration

V_lst = -70:1:40;
% V = @(t) -80 * (t < 10/1000)  - 40 * (t >= 10/1000) .* (t < 500/1000) - 80 * (t >= 500/1000);
t_induce = 1000;
J_Ca_lst_CaV12 = 0 * V_lst;
J_Ca_lst_CaV13 = 0 * V_lst;
I_NMDA_lst = 0 * V_lst;
parfor ind_V = 1:length(V_lst)
    V_induce = V_lst(ind_V);
    V = @(t) Vs * (t < t_induce/1000) + V_induce * (t >= t_induce/1000);
    
    t_end = 10;
    CaV12_test = @(t, y) CaV_model_Mahajan(y, V(t));
    sol_CaV12 = ode45(CaV12_test, [0, t_end], [c_Ca_SS;c_Ca;y0_CaV12]);
    t_eval = sol_CaV12.x; y_eval = sol_CaV12.y;
    Po = 1 - sum(y_eval(3:8, :), 1);
    g_CaV12 = 820; J_Ca = 0 * t_eval;
    for i = 1:length(t_eval)
        J_Ca(i) = CaV_flux(Po(i), g_CaV12, y_eval(2,i), V(t_eval(i)));
    end
    J_Ca_max_CaV12 = max(J_Ca); J_Ca_lst_CaV12(ind_V) = J_Ca_max_CaV12;
    
    CaV13_test = @(t, y) CaV13_model_Mahajan(y, V(t));
    sol_CaV13 = ode45(CaV13_test, [0, t_end], [c_Ca_SS;c_Ca;y0_CaV13]);
    t_eval = sol_CaV13.x; y_eval = sol_CaV13.y;
    Po = 1 - sum(y_eval(3:8, :), 1);
    g_CaV13 = 200; J_Ca = 0 * t_eval;
    for i = 1:length(t_eval)
        J_Ca(i) = CaV_flux(Po(i), g_CaV13, y_eval(2,i), V(t_eval(i)));
    end
    J_Ca_max_CaV13 = max(J_Ca); J_Ca_lst_CaV13(ind_V) = J_Ca_max_CaV13;
    
    NMDA_test = @(t, y) NMDA(s_NMDA, 1);
    sol_NMDA = ode45(NMDA_test, [0, t_end], [0]);
    t_eval = sol_NMDA.x; y_eval = sol_NMDA.y;
    I_NMDA_trial = 0 * t_eval;
    for i = 1:length(t_eval)
        [I_NMDA, I_NMDA_Ca] = NMDA_current(y_eval(i), V(t_eval(i)));
        I_NMDA_trial(i) = I_NMDA;
    end
    I_NMDA_max = max(I_NMDA_trial); I_NMDA_lst(ind_V) = I_NMDA_max;
end

figure();
plot(V_lst, -J_Ca_lst_CaV12/max(J_Ca_lst_CaV12));
hold on;
plot(V_lst, -J_Ca_lst_CaV13/max(J_Ca_lst_CaV13));
plot(V_lst(V_lst<0), -I_NMDA_lst(V_lst<0)/max(I_NMDA_lst));

legend('CaV1.2', 'CaV1.3','NMDAR');