% figure_counter = 4;
artifact = 'decay';
component = 8;

% spectrum
x = 1:lwdata.header.xstep:45;
y = squeeze(lwdata.data(1, component, 1, 1, 1, :));
y = log10(y);

fig = figure(figure_counter);
plot(x, y, 'linewidth', 3, 'color', [0.07 0.63 1])
xlim([-1, 47])
xlabel('frequency (Hz)')
ylabel('log10(\muV^2/Hz)')
set(gca, 'FontSize', 14)

filename = ['C:\Users\uzivatel\UCL\O365G-NOCIONS - dsulcova\GABA-AD_results\GABA_YC_figures\ICA\' artifact '_spectrum'];
savefig([filename '.fig'])
saveas(fig, [filename '.svg'], 'svg')
figure_counter = figure_counter + 1;

% timecourse
time_window = [-50, 300];
x = time_window(1):lwdata.header.xstep*1000:time_window(2);
x_start = (time_window(1)/1000 - lwdata.header.xstart)/lwdata.header.xstep;
x_end = (time_window(2)/1000 - lwdata.header.xstart)/lwdata.header.xstep;
y = squeeze(lwdata.data(1, component, 1, 1, 1, x_start:x_end));

fig = figure(figure_counter);
hold on

plot(x, y, 'linewidth', 3, 'color', [0.47 0.67 0.18])
% ylim([-60 15])
yl = get(gca, 'ylim');
ylim(yl)

rectangle('Position', [-5, yl(1), 15, yl(2) - yl(1)], 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none')
line([0, 0], yl, 'Color', [0 0 0], 'LineWidth', 3, 'LineStyle', '--')
child = get(gca,'Children');
set(gca,'Children',[child(3) child(1) child(2)])

xlim([-65, 315])
xlabel('time (ms)')
ylabel('amplitude (\muV)')
set(gca, 'FontSize', 14)
set(gca, 'Layer', 'Top')

filename = ['C:\Users\uzivatel\UCL\O365G-NOCIONS - dsulcova\GABA-AD_results\GABA_YC_figures\ICA\' artifact '_timecourse'];
savefig([filename '.fig'])
saveas(fig, [filename '.svg'], 'svg')
figure_counter = figure_counter + 1;
