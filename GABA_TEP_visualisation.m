%% parameters
clear all, clc

% dataset
target = 'AG';
group = 'YC';
participant = 1:20;
medication = {'placebo' 'alprazolam'}; 
time = {'pre' 'post'};
electrode = {'Cz'};
prefix = 'avg bl icfilt ica26 crop but fft-notchfilt prefilt prea28 visual reref ds art-sup bl ep dc chan-select'; 

% import a random header
load([prefix ' ' group num2str(participant(1)) ' ' medication{1} ' ' target ' pre.lw6'], '-mat')

% process
time_window = [-0.005, 0.3];
load('labels.mat'); labels = labels_TS(1:30); clear labels_TS labels_CS
z = 1.96;
x = [time_window(1):header.xstep:time_window(2)];
x_start = (time_window(1) - header.xstart)/header.xstep;
x_end = (time_window(2) - header.xstart)/header.xstep;
load('colours.mat'); col = colours([2 1 4 6],:);
alpha = 0.2;
figure_counter = 1;

%% AG TEP data extraction 
for m = 1:length(medication)
    for t = 1:length(time)
        for e = 1:length(electrode)
            for p = 1:length(participant)
                % load the dataset  
                load([prefix ' ' group num2str(participant(p)) ' ' medication{m} ' ' target ' ' time{t} '.mat']); 

                % choose the electrode                
                position = find(contains(labels, electrode{e}));
                AG_TEP(m, t, e, p, :) = squeeze(data(:, position, :, :, :, x_start:x_end))';
            end
        end
    end
end

% check matrix length
if size(AG_TEP, 5) ~= length(x)
    disp('ATTENTION: the length of x does not correspond to data.')
end

% prepare baseline statistics
for m = 1:length(medication)
    for t = 1:length(time)
        for e = 1:length(electrode)
            for i = 1:size(AG_TEP, 5)
                AG_TEP_mean(m, t, e, i) = mean(squeeze(AG_TEP(m, t, e, :, i)));
                AG_TEP_CI(m, t, e, i) = (std(squeeze(AG_TEP(m, t, e, :, i)))/sqrt(length(participant))) * z;
            end
        end
    end
end

%% plot baseline data 
% medication separately
data_visual = []; CI_visual = [];
for e = 1:length(electrode)
    % prepare data
    data_visual = squeeze(AG_TEP_mean(:, 1, e, :));
    CI_visual = squeeze(AG_TEP_CI(:, 1, e, :));

    % set limits of the figure
    fig = figure(figure_counter);
    hold on
    plot(x, data_visual(1, :) + CI_visual(1, :), 'b:', 'LineWidth', 0.5)
    plot(x, data_visual(1, :) - CI_visual(1, :), 'b:', 'LineWidth', 0.5)
    yl = get(gca, 'ylim');
    clf

    % launch the figure  
    p1 = plot(x, data_visual(1, :), 'b:', 'LineWidth', 0.5);
    rectangle('Position', [0, yl(1), 0.01, yl(2) - yl(1)], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')
    hold on

    % plot data
    for m = 1:length(medication)
        pl(m) = plot(x, data_visual(m, :), 'Color', col(m, :), 'LineWidth', 2.5);
        fi(m) = fill([x fliplr(x)],[data_visual(m, :) + CI_visual(m, :) fliplr(data_visual(m, :) - CI_visual(m, :))], ...
        col(m, :), 'FaceAlpha', alpha, 'linestyle', 'none');
        hold on
    end
    
    % add parameters
    xlim(time_window)
    ylim(yl)
    xlabel('time (s)')
    ylabel('amplitude (\muV)')
    set(gca, 'FontSize', 16)
    line([0, 0], get(gca,'ylim'), 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)
    lgd = legend(pl, medication, 'Location', 'southeast');
    lgd.FontSize = 16;
    hold  off
    
    % save figure    
    figure_name = ['TEP_AG_bl_' electrode{e}];
    savefig([figure_name])
    saveas(fig, [figure_name '.png'])
    pause(3)

    % update counter
    figure_counter = figure_counter + 1;
end

% both medications together
data_visual = []; CI_visual = [];
for e = 1:length(electrode)
    % prepare data
    for i = 1:size(AG_TEP, 4)
        data_visual(i) = mean(AG_TEP_bl(:, 1, e, i));          
        CI_visual(i) = mean(AG_TEP_CI(:, 1, e, i));
    end

    % set limits of the figure
    fig = figure(figure_counter);
    hold on
    plot(x, data_visual + CI_visual, 'b:', 'LineWidth', 0.5)
    plot(x, data_visual - CI_visual, 'b:', 'LineWidth', 0.5)
    yl = get(gca, 'ylim');
    clf

    % launch the figure  
    p1 = plot(x, data_visual, 'b:', 'LineWidth', 0.5);
    rectangle('Position', [0, yl(1), 0.01, yl(2) - yl(1)], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')
    hold on

    % plot data
    pl = plot(x, data_visual, 'Color', col(3, :), 'LineWidth', 2.5);
    fi = fill([x fliplr(x)],[data_visual + CI_visual fliplr(data_visual - CI_visual)], ...
    col(3, :), 'FaceAlpha', alpha, 'linestyle', 'none');
    hold on

    
    % add parameters
    xlim(time_window)
    ylim(yl)
    xlabel('time (s)')
    ylabel('amplitude (\muV)')
    set(gca, 'FontSize', 16)
    line([0, 0], get(gca,'ylim'), 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)
    hold  off
    
    % save figure    
    figure_name = ['TEP_AG_bl_both_' electrode{e}];
    savefig([figure_name])
    saveas(fig, [figure_name '.png'])
    pause(3)

    % update counter
    figure_counter = figure_counter + 1;
end

%% plot pre x post data 
% medication separately
data_visual = []; CI_visual = [];
for m = 1:length(medication)
    for e = 1:length(electrode)
        % prepare data
        data_visual = squeeze(AG_TEP_mean(m, :, e, :));
        CI_visual = squeeze(AG_TEP_CI(m, :, e, :));

        % set limits of the figure
        fig = figure(figure_counter);
        hold on
        plot(x, data_visual(1, :) + CI_visual(1, :), 'b:', 'LineWidth', 0.5)
        plot(x, data_visual(1, :) - CI_visual(1, :), 'b:', 'LineWidth', 0.5)
        yl = get(gca, 'ylim');
        clf

        % launch the figure  
        p1 = plot(x, data_visual(1, :), 'b:', 'LineWidth', 0.5);
        rectangle('Position', [0, yl(1), 0.01, yl(2) - yl(1)], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')
        hold on

        % plot data
        for t = 1:length(time)
            pl(t) = plot(x, data_visual(t, :), 'Color', col((m-1) * 2 + t, :), 'LineWidth', 2.5);
            fi(t) = fill([x fliplr(x)],[data_visual(t, :) + CI_visual(t, :) fliplr(data_visual(t, :) - CI_visual(t, :))], ...
            col((m-1) * 2 + t, :), 'FaceAlpha', alpha, 'linestyle', 'none');
            hold on
        end

        % add parameters
        xlim(time_window)
        ylim(yl)
        xlabel('time (s)')
        ylabel('amplitude (\muV)')
        set(gca, 'FontSize', 16)
        line([0, 0], get(gca,'ylim'), 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)
        lgd = legend(pl, {'pre medication' 'post medication'}, 'Location', 'southeast');
        lgd.FontSize = 16;
        hold  off

        % save figure    
        figure_name = ['TEP_' target '_' medication{m} '_' electrode{e}];
        savefig([figure_name])
        saveas(fig, [figure_name '.png'])
        pause(3)

        % update counter
        figure_counter = figure_counter + 1;
    end
end