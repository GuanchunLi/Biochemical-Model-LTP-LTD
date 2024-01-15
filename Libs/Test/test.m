freq_cal = 10; t_induce = 10; duration = 2;
cal_source = @(t) 0.05 + 4 * sum(exp(-(t-t_induce-(1:(min(t-t_induce, duration)*freq_cal))/freq_cal)/0.1));
freq_cal = 60;
cal_source_2 = @(t) 0.05 + 4 * sum(exp(-(t-t_induce-(1:(min(t-t_induce, duration)*freq_cal))/freq_cal)/0.1));

t_lst = 5:1e-3:20;
cal_lst = zeros(size(t_lst));
cal_lst_2 = zeros(size(t_lst));
for i = 1:length(cal_lst)
    cal_lst(i) = cal_source(t_lst(i));
    cal_lst_2(i) = cal_source_2(t_lst(i));
end

fig = figure(1);
plot(t_lst, cal_lst); 

hold on;
plot(t_lst, cal_lst_2); 

xlabel('Time (s)');
ylabel('Calcium (\mu M)');
legend('10Hz, 0.05 \muM, 2s', '60Hz, 0.5 \muM, 2s', 'Location', 'east');
% saveas(fig, ['Result/calcium_dynamics_two.png']);

% title(['Calcium dynamics when stimulus = ', num2str(freq_cal), 'Hz']);
% saveas(fig, ['../Result/Ca_dynamics_', num2str(freq_cal), 'Hz.png']);