%% 
figh = figure;
tiledlayout(figh, 4, 6);

nexttile(1, [1 3])
plot(t_spline, V_Xe_reg(:,:,1), 'LineWidth', 2)
xlabel('Time [s]','Interpreter','latex')
ylabel('Force [N]','Interpreter','latex')
hold on
yline(Par_PDF.thres_meas, "r-", "threshold", "LineWidth", 2.5);
hold off
grid on
legend("$V(t)$")
xlim([0, 10])

set(findall(figh,'-property','Interpreter'),'Interpreter','latex') 
set(findall(figh,'-property','TickLabelInterpreter'), ...
    'TickLabelInterpreter','latex')
set(findall(figh, "-property", "FontSize"), "FontSize", 18)