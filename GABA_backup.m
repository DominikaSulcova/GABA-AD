%% 5) DISS
% m = 1; t = 1; s = 1, p = 1; i = 1;
% normalize data by GMFP
for m = 1:length(medication)
    for t = 1:length(time)
        for s = 1:length(stimulus) 
            for p = 1:length(participant)
                for i = 1:size(GABA_data, 6)
                    % devide data at each time point by GMFP
                    GABA_data_norm(m, t, s, p, :, i) = squeeze(GABA_data(m, t, s, p, 1:30, i)) / GABA_GMFP_subject(m, t, s, p, i);
                end
            end
        end
    end
end
clear m t s p i

% ----- baseline -----
% select baseline data of both sessions
for s = 1:length(stimulus)
    counter = 1;
    for p = 1:length(participant)
        for m = 1:length(medication)
            data_temp(s, counter, :, :) = squeeze(GABA_data_norm(m, 1, s, p, :, :));
            counter = counter + 1;
        end
    end
end
clear s p m counter

% average across subjects
for s = 1:length(stimulus)
    data(s, :, :) = squeeze(mean(data_temp(s, :, :, :), 2));
end
clear s data_temp

% calculate DISS and spatial correlation C
GABA_DISS = struct;
for i = 1:size(GABA_data_norm, 6)
    % between TS and CS stimuli
    diff = squeeze(data(2, :, i) - data(1, :, i));
    GABA_DISS.baseline.stimulus(1, i) = sqrt(mean(diff.^2));
    GABA_DISS.baseline.stimulus(2, i) = 1 - (GABA_DISS.baseline.stimulus(1, i)^2)/2;

    % between ppTMS and TS stimuli
    diff = squeeze(data(3, :, i) - data(2, :, i));
    GABA_DISS.baseline.SICI(1, i) = sqrt(mean(diff.^2));
    GABA_DISS.baseline.SICI(2, i) = 1 - (GABA_DISS.baseline.SICI(1, i)^2)/2;
end
clear p i diff

% plot baseline TS-CS
data_visual = GABA_DISS.baseline.stimulus;
fig = figure(figure_counter);
hold on
title('Global dissimilarity: baseline TS - CS', 'FontSize', 16, 'FontWeight', 'bold')
set(gca, 'fontsize', 12); xlabel('time (s)')
ax = gca; ax.XColor = [0.5 0.5 0.5];
plot(x, data_visual(1, :), 'Color', colours(2, :), 'LineWidth', 2.5)
fill([x fliplr(x)],[data_visual(2, :) zeros(1, length(data_visual(2, :)))], ...
    colours(1, :) , 'facealpha', 0.5, 'linestyle', 'none');
line(time_window, [0 0], 'Color', [0, 0, 0], 'LineWidth', 0.5)
line(time_window, [1 1], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5)
line([0, 0], [-0.75 1.75], 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5) 
ylim([-0.75 1.75]); 
lgd = legend({'DISS index' 'spatial correlation'}, 'Location', 'southeast');
lgd.FontSize = 12;
hold off
figure_name = ['DISS_' target '_bl_stimulus'];
savefig([folderpath '\' filename '_figures\' figure_name '.fig'])
saveas(fig, [folderpath '\' filename '_figures\' figure_name '.png'])
figure_counter = figure_counter + 1;
clear data_visual fig ax lgd figure_name

% plot baseline ppTMS-TS (SICI)
data_visual = GABA_DISS.baseline.SICI;
fig = figure(figure_counter);
hold on
title('Global dissimilarity: baseline SICI', 'FontSize', 16, 'FontWeight', 'bold')
set(gca, 'fontsize', 12); xlabel('time (s)')
ax = gca; ax.XColor = [0.5 0.5 0.5];
plot(x, data_visual(1, :), 'Color', colours(6, :), 'LineWidth', 2.5)
fill([x fliplr(x)],[data_visual(2, :) zeros(1, length(data_visual(2, :)))], ...
    colours(5, :) , 'facealpha', 0.5, 'linestyle', 'none');
line(time_window, [0 0], 'Color', [0, 0, 0], 'LineWidth', 0.5)
line(time_window, [1 1], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5)
line([0, 0], [-0.75 1.75], 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5) 
ylim([-0.75 1.75]); 
lgd = legend({'DISS index' 'spatial correlation'}, 'Location', 'southeast');
lgd.FontSize = 12;
hold off
figure_name = ['DISS_' target '_bl_SICI'];
savefig([folderpath '\' filename '_figures\' figure_name '.fig'])
saveas(fig, [folderpath '\' filename '_figures\' figure_name '.png'])
figure_counter = figure_counter + 1;
clear data_visual fig ax lgd figure_name

% ----- pre-post medication -----
% calculate DISS and spatial correlation C
for m = 1:length(medication)
    for s = 1:length(stimulus)
        for i = 1:size(GABA_data_norm, 6)
            diff = mean(squeeze(GABA_data_norm(m, 2, s, :, :, i)), 1) - mean(squeeze(GABA_data_norm(m, 1, s, :, :, i)), 1);
            GABA_DISS.medication(1, m, s, i) = sqrt(mean(diff.^2));
            GABA_DISS.medication(2, m, s, i) = 1 - (GABA_DISS.medication(1, m, s, i)^2)/2;
        end
    end
end
clear m s i diff

% plot 
for m = 1:length(medication)
    for s = 1:length(stimulus) 
        data_visual = squeeze(GABA_DISS.medication(:, m, s, :));
        fig = figure(figure_counter);
        hold on
        title(sprintf('Global dissimilarity: %s, %s', medication{m}, stimulus{s}), 'FontSize', 16, 'FontWeight', 'bold')
        set(gca, 'fontsize', 12); xlabel('time (s)')
        ax = gca; ax.XColor = [0.5 0.5 0.5];
        plot(x, data_visual(1, :), 'Color', colours((m-1)*2 + 2, :), 'LineWidth', 2.5)
        fill([x fliplr(x)],[data_visual(2, :) zeros(1, length(data_visual(2, :)))], ...
            colours((m-1)*2 + 1, :) , 'facealpha', 0.5, 'linestyle', 'none');
        line(time_window, [0 0], 'Color', [0, 0, 0], 'LineWidth', 0.5)
        line(time_window, [1 1], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5)
        line([0, 0], [-0.75 1.75], 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5) 
        ylim([-0.75 1.75]); 
        lgd = legend({'DISS index' 'spatial correlation'}, 'Location', 'southeast');
        lgd.FontSize = 12;
        hold off
        figure_name = ['DISS_' target '_' medication{m} '_' stimulus{s}];
        savefig([folderpath '\' filename '_figures\' figure_name '.fig'])
        saveas(fig, [folderpath '\' filename '_figures\' figure_name '.png'])
        figure_counter = figure_counter + 1;
    end
end
clear data_visual fig ax lgd figure_name