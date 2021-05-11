clear all; clc;

% load data
load('TEPs_final.mat');
TEPs = TEPs_final;

% parameters
participant = 1:20;
peaks = {'P30' 'N45' 'P60' 'N100' 'P180'};
medication = {'placebo' 'alprazolam'};
time = {'pre' 'post'};
stimulus = {'CS'};
load('colours.mat');

%% sp-TMS TEP visualization
figure_counter = 1;

for k = 1:length(peaks)
    % extract data
    rows = (categorical(TEPs.peak) == peaks{k} & categorical(TEPs.stimulus) == char(stimulus{1}));
    data = TEPs(rows, :);
    height(data)

    %% mean TEP amplitude: absolute values
    % calculate group average
    for a = 1:length(medication)
        for b = 1:length(time)
            avg_amp(a, b) = mean(data.amplitude(strcmp(data.medication, medication{a}) & strcmp(data.time, time{b})));
            avg_amp_std(a, b) = std(data.amplitude(strcmp(data.medication, medication{a}) & strcmp(data.time, time{b})));
            avg_amp_sem(a, b) = avg_amp_std(a, b)/sqrt(length(participant));
        end
    end

    % plot the bar graph with both medications 
    fig = figure(figure_counter);
    bar(avg_amp', 'EdgeColor', 'none')
    colormap(colours)
    hold on

    % add errorbars
    ngroups = length(time);
    nbars = length(medication);
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, avg_amp(i, :), avg_amp_sem(i, :), 'k', 'linestyle', 'none');
    end

    % set figure names 
    fig_name = [char(stimulus{1}) ' - ' data.peak{1} ': mean amplitude'];
    peak_name = data.peak{1};
    
    % set other parameters
    title(fig_name, 'FontWeight', 'bold', 'FontSize', 20)
    legend('placebo', 'alprazolam', 'Location', 'bestoutside')
    xlabel('time to medication')
    ylabel('amplitude (?V)')
    set(gca, 'xtick', 1:length(time), 'xticklabel', time)
    set(gca, 'Fontsize', 16)
    hold on

    % save and close figure
    figure_counter = figure_counter + 1;
    savefig([char(stimulus{1}) '_' peak_name '_absolute.fig'])
    pause(5)
    close(fig)
    
%     %% mean TEP amplitude: normalized values
%     % calculate normalized amplitudes
%     for p = 1:length(participant)
%         % placebo
%         rows_pre = (categorical(data.medication) == 'placebo' & data.subject == participant(p) & categorical(data.time) == 'pre');
%         rows_post = (categorical(data.medication) == 'placebo' & data.subject == participant(p) & categorical(data.time) == 'post');
%         placebo_change(p) = data.amplitude(rows_post)/data.amplitude(rows_pre) * 100; 
%         
%         % alprazolam
%         rows_pre = (categorical(data.medication) == 'alprazolam' & data.subject == participant(p) & categorical(data.time) == 'pre');
%         rows_post = (categorical(data.medication) == 'alprazolam' & data.subject == participant(p) & categorical(data.time) == 'post');
%         alprazolam_change(p) = data.amplitude(rows_post)/data.amplitude(rows_pre) * 100;         
%     end
%     
%     % group average amplitude
%     avg_amp(1:2, 1) = 100;
%     avg_amp(1, 2) = mean(placebo_change);
%     avg_amp(2, 2) = mean(alprazolam_change);
%     
%     % group average sem
%     avg_amp_sem(1:2, 1) = 0;
%     avg_amp_sem(1, 2) = std(placebo_change)/sqrt(length(participant));
%     avg_amp_sem(2, 2) = std(alprazolam_change)/sqrt(length(participant));
%     
%     % plot the bar graph with both medications 
%     fig = figure(figure_counter);
%     bar(avg_amp', 'EdgeColor', 'none')
%     colormap(colours)
%     hold on
% 
%     % add errorbars
%     ngroups = length(time);
%     nbars = length(medication);
%     groupwidth = min(0.8, nbars/(nbars + 1.5));
%     for i = 1:nbars
%         x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%         errorbar(x, avg_amp(i, :), avg_amp_sem(i, :), 'k', 'linestyle', 'none');
%     end
%     
%     % add baseline line
%     x_pos = get(gca,'xlim');
%     plot(x_pos, [100 100], 'LineStyle', ':', 'Color', [0 0 0]);
%     hold on
% 
%     % set figure names 
%     if k == 2
%         fig_name = [data.peak{1} ' - central EOIs: normalized mean amplitude'];
%         peak_name = [data.peak{1} '_central'];
%     elseif k == 3
%         fig_name = [data.peak{1} ' - lateral EOIs: normalized mean amplitude'];
%         peak_name = [data.peak{1} '_lateral'];
%     else
%         fig_name = [data.peak{1} ': normalized mean amplitude'];
%         peak_name = data.peak{1};
%     end
%     
%     % set other parameters
%     title(fig_name, 'FontWeight', 'bold', 'FontSize', 20)
%     legend('placebo', 'alprazolam', 'Location', 'bestoutside')
%     xlabel('time to medication')
%     ylabel('amplitude (% baseline)')
%     set(gca, 'xtick', 1:length(time), 'xticklabel', time)
%     set(gca, 'Fontsize', 16)
%     hold on
% 
%     % save and close figure
%     figure_counter = figure_counter + 1;
%     savefig([peak_name '_normalized.fig'])
%     close(fig)
end

%% baseline SICI - visualization
for k = 1:length(electrode)
    % extract data
    rows = (categorical(TEPs.channel) == electrode{k} & categorical(TEPs.time) == 'pre');
    data = TEPs(rows, :);
    height(data)

    %% mean TEP amplitude: absolute values
    % calculate group average
    for a = 1:length(medication)
        for b = 1:length(stimulus)
            avg_amp(a, b) = mean(data.amplitude(strcmp(data.medication, medication{a}) & strcmp(data.stimulus, stimulus{b})));
            avg_amp_std(a, b) = std(data.amplitude(strcmp(data.medication, medication{a}) & strcmp(data.stimulus, stimulus{b})));
            avg_amp_sem(a, b) = avg_amp_std(a, b)/sqrt(length(participant));
        end
    end

    % plot the bar graph with both medications 
    fig = figure(figure_counter);
    bar(avg_amp', 'EdgeColor', 'none')
    colormap(colours)
    hold on

    % add errorbars
    ngroups = length(stimulus);
    nbars = length(medication);
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, avg_amp(i, :), avg_amp_sem(i, :), 'k', 'linestyle', 'none');
    end

    % set figure names 
    if k == 2
        fig_name = ['SICI - baseline:' data.peak{1} ' - central EOIs: absolute amplitude'];
        peak_name = [data.peak{1} '_central'];
    elseif k == 3
        fig_name = ['SICI - baseline:' data.peak{1} ' - lateral EOIs: absolute amplitude'];
        peak_name = [data.peak{1} '_lateral'];
    else
        fig_name = ['SICI - baseline:' data.peak{1} ': absolute amplitude'];
        peak_name = data.peak{1};
    end
    
    % set other parameters
    title(fig_name, 'FontWeight', 'bold', 'FontSize', 20)
    legend('placebo', 'alprazolam', 'Location', 'bestoutside')
    xlabel('stimulus')
    ylabel('amplitude (?V)')
    set(gca, 'xtick', 1:length(time), 'xticklabel', stimulus)
    set(gca, 'Fontsize', 16)
    hold on

    % save and close figure
    figure_counter = figure_counter + 1;
    savefig(['SICI_baseline_' peak_name '_absolute.fig'])
    close(fig)
    
    %% mean TEP amplitude: normalized values
    % calculate normalized amplitudes
    for p = 1:length(participant)
        % placebo
        rows_TS = (categorical(data.medication) == 'placebo' & data.subject == participant(p) & categorical(data.stimulus) == 'TS');
        rows_ppTMS = (categorical(data.medication) == 'placebo' & data.subject == participant(p) & categorical(data.stimulus) == 'ppTMS');
        placebo_change(p) = data.amplitude(rows_ppTMS)/data.amplitude(rows_TS) * 100; 
        
        % alprazolam
        rows_TS = (categorical(data.medication) == 'alprazolam' & data.subject == participant(p) & categorical(data.stimulus) == 'TS');
        rows_ppTMS = (categorical(data.medication) == 'alprazolam' & data.subject == participant(p) & categorical(data.stimulus) == 'ppTMS');
        alprazolam_change(p) = data.amplitude(rows_ppTMS)/data.amplitude(rows_TS) * 100;         
    end
    
    % group average amplitude
    avg_amp(1:2, 1) = 100;
    avg_amp(1, 2) = mean(placebo_change);
    avg_amp(2, 2) = mean(alprazolam_change);
    
    % group average sem
    avg_amp_sem(1:2, 1) = 0;
    avg_amp_sem(1, 2) = std(placebo_change)/sqrt(length(participant));
    avg_amp_sem(2, 2) = std(alprazolam_change)/sqrt(length(participant));
    
    % plot the bar graph with both medications 
    fig = figure(figure_counter);
    bar(avg_amp', 'EdgeColor', 'none')
    colormap(colours)
    hold on

    % add errorbars
    ngroups = length(stimulus);
    nbars = length(medication);
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, avg_amp(i, :), avg_amp_sem(i, :), 'k', 'linestyle', 'none');
    end
    
    % add baseline line
    x_pos = get(gca,'xlim');
    plot(x_pos, [100 100], 'LineStyle', ':', 'Color', [0 0 0]);
    hold on

    % set figure names 
    if k == 2
        fig_name = ['SICI - baseline:' data.peak{1} ' - central EOIs: normalized amplitude'];
        peak_name = [data.peak{1} '_central'];
    elseif k == 3
        fig_name = ['SICI - baseline:' data.peak{1} ' - lateral EOIs: normalized amplitude'];
        peak_name = [data.peak{1} '_lateral'];
    else
        fig_name = ['SICI - baseline:' data.peak{1} ': normalized amplitude'];
        peak_name = data.peak{1};
    end
    
    % set other parameters
    title(fig_name, 'FontWeight', 'bold', 'FontSize', 20)
    legend('placebo', 'alprazolam', 'Location', 'bestoutside')
    xlabel('stimulus')
    ylabel('amplitude (% baseline)')
    set(gca, 'xtick', 1:length(stimulus), 'xticklabel', stimulus)
    set(gca, 'Fontsize', 16)
    hold on

    % save and close figure
    figure_counter = figure_counter + 1;
    savefig(['SICI_baseline_' peak_name '_normalized.fig'])
    close(fig)
end

