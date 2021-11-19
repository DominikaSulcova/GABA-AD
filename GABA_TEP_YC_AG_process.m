%% GABA-AD: TMS-EVOKED POTENTIALS - ANGULAR GYRUS
% Written by Dominika for GABA-AD project (2020-21)

%% parameters
clear all; clc

% ----- adjustable parameters -----
% dataset
prefix = 'GABA';
group = 'YC';
target = 'AG';
participant = 1:20;
medication = {'placebo' 'alprazolam'}; 
time = {'pre' 'post'};

% visualization 
time_window = [-0.05, 0.3];
z = 1.96;
alpha = 0.2;
% --------------------------------

% navigate to the Git folder --> source of saved default files
folder_git = uigetdir(pwd, 'Choose the Git folder');

% load default header
load([folder_git '\GABA_header_default.mat'])
for l = 1:size(header.chanlocs, 2)
    labels{l} = header.chanlocs(l).labels;
end
clear l

% visualization calculated params
figure_counter = 1;
xstep = header.xstep; 
xstart = header.xstart;
x = [time_window(1):xstep:time_window(2)];

% check for colour scheme
answer = questdlg('Do you want to choose a new colour scheme?', 'Colour scheme', 'YES', 'NO', 'NO'); 
switch answer
    case 'YES'
        a = 1;
        for p = 1:6
            for c = 1:length(current)
               colours(a, :) = uisetcolor; 
               a = a + 1;
            end
        end
        save('colours.mat', 'colours'); 
    case 'NO'
    load([folder_git '\GABA_colours.mat'])
end
clear a answer

% create output folders
filename = ['GABA_' group '_' target '_TEPs'];
folder_results = uigetdir(pwd, 'Choose the Results folder');
output_file = [folder_results '\' filename '.mat'];
folder_export = [folder_results '\GABA_' group '_export'];
if ~exist(folder_export)
    mkdir(folder_export)
end
folder_figures = [folder_results '\GABA_' group '_figures'];
if ~exist(folder_figures)
    mkdir(folder_figures)  
end
clear folder_results

%% 1) extract individual data
% load data 
for m = 1:length(medication)
    for t = 1:length(time)
        for p = 1:length(participant)
           % define participant
            if p < 10
                subj = ['0' num2str(participant(p))];
            else
                subj = num2str(participant(p));
            end

            % identify dataset
            dataset_name = [prefix ' ' group ' ' subj ' ' target ' ' medication{m} ' ' time{t}];

            % identify cropping limits
            if m == 1 && t == 1 && p == 1
                % calculate limits
                load([dataset_name '.lw6'], '-mat')                    
                x_start = (time_window(1) - header.xstart)/xstep;
                x_end = (time_window(2) - header.xstart)/xstep;

                % reload the default header
                load([folder_git '\GABA_header_default.mat'])
            end

            % load individual dataset
            load([dataset_name '.mat'])

            % append the data in the data matrix
            GABA_data(m, t, p, :, :) = squeeze(data(:, :, :, :, :, x_start:x_end));                                                        
        end
    end
end
clear m t p data subj
disp(['Data import finished. Datasize: ' num2str(size(GABA_data))])

% save dataset to the global MATLAB file
if ~exist(output_file) 
    save(output_file, 'GABA_data');
else
    save(output_file, 'GABA_data', '-append');
end

%% 2) preliminary visualization 
% ----- decide output parameters -----
electrode = {'Cz'};
% ------------------------------------
% load data if not loaded
if exist('GABA_data') ~= 1
    load([folder_results '\' filename '.mat'], 'GABA_data')
end

% prepare average statistics
for m = 1:length(medication)
    for t = 1:length(time)
        for e = 1:size(GABA_data, 4)
            for i = 1:size(GABA_data, 5)
                GABA_data_mean(m, t, e, i) = mean(squeeze(GABA_data(m, t, :, e, i)));
                GABA_data_CI(m, t, e, i) = (std(squeeze(GABA_data(m, t, :, e, i)))/sqrt(length(participant))) * z;
            end
        end
    end
end
clear m t e i

% plot the baseline from both sessions
for e = 1:length(electrode) 
    % identify the electrode
    e_n = find(contains(labels, electrode{e}));

    % prepare data
    data_visual = squeeze(GABA_data_mean(:, 1, e_n, :));
    CI_visual = squeeze(GABA_data_CI(:, 1, e_n, :));

    % launch the figure
    fig = figure(figure_counter);
    hold on

    % set limits of the figure
    plot(x, data_visual(1, :) + CI_visual(1, :), 'b:', 'LineWidth', 0.5)
    plot(x, data_visual(1, :) - CI_visual(1, :), 'b:', 'LineWidth', 0.5)
    yl = get(gca, 'ylim'); yl(1) = yl(1) - 0.2; yl(2) = yl(2) + 0.3;
    clf, hold on

    % shade interpolated interval 
    plot(x, data_visual(1, :), 'b:', 'LineWidth', 0.5);
    rectangle('Position', [0, yl(1), 0.01, yl(2) - yl(1)], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')

    % loop through datasets to plot
    for m = 1:length(medication)        
        P(m) = plot(x, data_visual(m, :), 'Color', colours((m - 1)*2 + 1, :), 'LineWidth', 2.5);
        F(m) = fill([x fliplr(x)],[data_visual(m, :) + CI_visual(m, :) fliplr(data_visual(m, :) - CI_visual(m, :))], ...
            colours((m - 1)*2 + 1, :), 'FaceAlpha', alpha, 'linestyle', 'none');
    end

    % add other parameters
    title([target ' - ' electrode{e} ' electrode: baseline TEP'])
    xlabel('time (s)')
    ylabel('amplitude (\muV)')
    set(gca, 'FontSize', 14)
    xlim(time_window)
    ylim(yl)
    line([0, 0], get(gca,'ylim'), 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)
    lgd = legend(P, medication, 'Location', 'southeast');
    lgd.FontSize = 14;
    hold off

    % name and save figure
    figure_name = ['TEP_' target '_bl_' electrode{e}];
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.png'])

    % update figure counter
    figure_counter = figure_counter + 1 ;
end
clear m e e_n fig lgd data_visual CI_visual figure_name P F yl

% plot baseline vs. post-medication 
for e = 1:length(electrode) 
    for m = 1:length(medication)
        % identify the electrode
        e_n = find(contains(labels, electrode{e}));

        % prepare data
        data_visual = squeeze(GABA_data_mean(m, :, e_n, :));
        CI_visual = squeeze(GABA_data_CI(m, :, e_n, :));

        % launch the figure
        fig = figure(figure_counter);
        hold on

        % set limits of the figure
        plot(x, data_visual(1, :) + CI_visual(1, :), 'b:', 'LineWidth', 0.5)
        plot(x, data_visual(1, :) - CI_visual(1, :), 'b:', 'LineWidth', 0.5)
        yl = get(gca, 'ylim'); yl(1) = yl(1) - 0.2; yl(2) = yl(2) + 0.3;
        clf, hold on

        % shade interpolated interval 
        plot(x, data_visual(1, :), 'b:', 'LineWidth', 0.5);
        rectangle('Position', [0, yl(1), 0.01, yl(2) - yl(1)], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')

        % loop through datasets to plot
        for t = 1:length(time)        
            P(t) = plot(x, data_visual(t, :), 'Color', colours((m - 1)*2 + t, :), 'LineWidth', 2.5);
            F(t) = fill([x fliplr(x)],[data_visual(t, :) + CI_visual(t, :) fliplr(data_visual(t, :) - CI_visual(t, :))], ...
                colours((m - 1)*2 + t, :), 'FaceAlpha', alpha, 'linestyle', 'none');
        end

        % add other parameters
        title([target ' - ' electrode{e} ' electrode: ' medication{m}])
        xlabel('time (s)')
        ylabel('amplitude (\muV)')
        set(gca, 'FontSize', 14)
        xlim(time_window)
        ylim(yl)
        line([0, 0], get(gca,'ylim'), 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)
        lgd = legend(P, {'baseline' 'post medication'}, 'Location', 'southeast');
        lgd.FontSize = 14;
        hold off

        % name and save figure
        figure_name = ['TEP_' target '_' medication{m} '_' electrode{e}];
        savefig([folder_figures '\' figure_name '.fig'])
        saveas(fig, [folder_figures '\' figure_name '.png'])

        % update figure counter
        figure_counter = figure_counter + 1 ;
    end
end
clear m t e e_n fig lgd data_visual CI_visual figure_name P F yl

% plot baseline butterfly plots per session/condition
for m = 1:length(medication)
     % prepare data
    data_visual = squeeze(GABA_data_mean(m, 1, 1:30, :));

    % launch the figure
    fig = figure(figure_counter);
    hold on

    % set limits of the figure
    yl = [-4.5, 4.5];

    % shade interpolated interval
    plot(x, data_visual(1, :), 'b:', 'LineWidth', 0.5);
    rectangle('Position', [-0.005, yl(1), 0.015, yl(2) - yl(1)], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')

    % plot all channels
    for e = 1:size(data_visual, 1)
        plot(x, data_visual(e, :), 'Color', colours(4 + m, :))
    end

    % mark TMS stimulus and zerol line
    line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2)

    % add other parameters
    title(sprintf('%s - TEP: %s', target, medication{m}))
    ylabel('amplitude (\muV)')
    set(gca, 'FontSize', 14)
    xlim(time_window)
    ylim(yl)
    hold off

    % change figure size
    fig.Position = [250 250 500 350];

    % name and save figure
    figure_name = ['TEP_' target '_bl_' medication{m}];
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.png'])

    % update figure counter
    figure_counter = figure_counter + 1 ;
end
clear m e data_visual fig yl  

% save dataset to the global MATLAB file
save(output_file, 'GABA_data_mean', 'GABA_data_CI', '-append');
clear electrode 

%% 3) GFP 
% compute individual GFP
for m = 1:length(medication)
    for t = 1:length(time)
        for p = 1:length(participant)
            % calculate GFP (exclude target channel)
            GABA_GFP(m, t, p, :) = std(squeeze(GABA_data(m, t, p, 1:30, :)), 1);  
        end
    end
end
clear m t s p

% append new variable to the general MATLAB file
save(output_file, 'GABA_GFP', '-append');

%% 4) GFP - peak identification
% ----- decide output parameters -----
labeled = 'off';
max_peaks = 6;
% ------------------------------------
if exist('GABA_data') ~= 1
    load([folder_results '\' filename '.mat'], 'GABA_data')
end

% calculate mean GFP (exclude target channel)
for m = 1:length(medication)
    for t = 1:length(time)
        GABA_GFP_mean(m, t, :) = std(squeeze(GABA_data_mean(m, t, 1:30, :)), 1);  
    end 
end
clear m t 

% ----- plot global GFP -----
% pool all conditions together
for i = 1:size(GABA_GFP_mean, 3)
    data_visual(i) = mean(GABA_GFP_mean(:, :, i), 'all');
end
clear i

% launch the figure
fig = figure(figure_counter);
hold on

% extract peak latencies
h_axis(1) = subplot(3, max_peaks, [1 : 2*max_peaks]);
GABA_peaks = gfp_plot(x, data_visual, time_window, xstep, labeled, 'max_peaks', max_peaks);
title([target ': GFP - all conditions'], 'fontsize', 16, 'fontweight', 'bold')

% choose data for topoplots 
for e = 1:30
    for i = 1:size(GABA_data_mean, 4)
        data_topoplot(1, e, 1, 1, 1, i) = mean(GABA_data_mean(:, :, e, i), 'all');
    end
end
clear e

% add topoplots
for t = 1:length(GABA_peaks)
    % plot the topoplot
    h_axis(1 + t) = subplot(3, max_peaks, 2*max_peaks + t);
    topo_plot(header, data_topoplot, GABA_peaks(t), time_window(1), [-2, 2])

    % shift down
    pos = get(h_axis(1 + t), 'Position');
    pos(2) = pos(2) - 0.05;
    set(h_axis(1 + t), 'Position', pos);

    % add timing
    text(-0.3, -0.8, sprintf('%1.0f ms', GABA_peaks(t)*1000), 'Color', [1 0 0], 'FontSize', 14)
end
clear t
hold off

% save figure
figure_name = ['GFP_' target];
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.png'])

% update figure counteer
figure_counter = figure_counter + 1;
clear data_visual data_topoplot fig figure_name pos h_axis GABA_peaks

% append new variables to the general MATLAB file
save(output_file, 'GABA_GFP_mean', '-append');
clear data_visual data_topoplot fig figure_name pos h_axis GABA_peaks labeled max_peaks 

%% 5) GFP - plot difference
% ----- decide output parameters -----
TOI_peaks = [0.025 0.045 0.075 0.110 0.180];
peaks = {'P25' 'N45' 'P75' 'N100' 'P180'};
% ------------------------------------
% plot baseline vs. post-medication
lines = {':' '-'};
for m = 1:length(medication)
    % load grand average data                
    data_visual = squeeze(mean(squeeze(GABA_GFP(m, :, :, :)), 2));

    % launch the figure
    fig = figure(figure_counter);
    h_axis(1) = subplot(4, length(TOI_peaks) + 1 , 1 : 2*(length(TOI_peaks) + 1));
    title(sprintf('%s: GFP - %s', target, medication{m}),...
        'FontSize', 16, 'FontWeight', 'bold')
    set(gca, 'fontsize', 12)
    xlabel('time (s)')
    ylabel('GFP (\muV)')
    hold on

    % set limits of the figure
    plot(x, data_visual(1, :))
    yl = get(gca, 'ylim'); yl(1) = yl(1) - 0.2; yl(2) = yl(2) + 0.3;
    xl = get(gca, 'xlim');
    cla(h_axis(1))
    ylim(yl); xlim(xl);

    % plot interpolated part
    rectangle('Position', [0, yl(1)+0.01, 0.01, yl(2) - yl(1)],...
        'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none'); 

    % plot the shading
    F = fill([x fliplr(x)],[data_visual(1, :) fliplr(data_visual(2, :))], ...
        colours((m-1)*2 + 1, :), 'linestyle', 'none', 'facealpha', 0.5);
    F.Annotation.LegendInformation.IconDisplayStyle = 'off';

    % plot GFP
    for t = 1:length(time)
        % plot signal
        plot(x, data_visual(t, :), 'Parent', h_axis(1), ...
            'Color', colours((m-1)*2 + 2, :), 'LineStyle', lines{t}, 'LineWidth', 2.5);
    end

    % mark TMS stimulus
    line([0, 0], yl, 'Parent', h_axis(1), 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)  

    % add legend
    lgd = legend({'baseline' 'post medication'}, 'Location', 'southeast');
    lgd.FontSize = 12; 

    % add topoplots
    for t = 1:length(time)
        % choose data for topoplots 
        for e = 1:30
            for i = 1:size(GABA_data, 5)
                data_topoplot(1, e, 1, 1, 1, i) = mean(squeeze(GABA_data(m, t, :, e, i)), 'all');
            end
        end

        % add topoplots for all peaks
        for k = 1:length(TOI_peaks)
            n = 1 + (t-1)*length(TOI_peaks) + k;

            % plot the topoplot
            h_axis(n) = subplot(4, length(TOI_peaks) + 1, 2*(length(TOI_peaks) + 1) + (t-1)*(length(TOI_peaks)+1) + k);
            topo_plot(header, data_topoplot, TOI_peaks(k), time_window(1), [-2, 2])

            % modifies the layout
            if t == 1
                pos = get(h_axis(n), 'Position');
                pos(2) = pos(2) - 0.04;
                set(h_axis(n), 'Position', pos);
            else
                text(-0.2, -0.7, peaks{k}, 'FontWeight', 'bold', 'FontSize', 12)
            end

            if k == length(TOI_peaks)
                if t == 1
                    descript = 'baseline';
                else
                    descript = sprintf('post \nmedication');
                end                    
                text(0.8, 0, descript, 'Parent', h_axis(n), ...
                    'Color', [0 0 0], 'FontSize', 12)
            end
        end
    end
    hold off

    % change figure position
    pos = get(gcf, 'Position');
    fig.Position = [pos(1) pos(2)-200 660 550];

    % save figure
    figure_name = ['GFP_' target '_' medication{m}];
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.png'])

    % update figure counteer
    figure_counter = figure_counter + 1;
end
clear peaks lines m t fig h_axis yl xl F lgd pos k n  P R L descript figure_name e i

%% 6) DISS
% load data if not loaded
if exist('GABA_GFP') ~= 1
    load([folder_results '\' filename '.mat'], 'GABA_GFP_subject')
end

% normalize data by GFP
for m = 1:length(medication)
    for t = 1:length(time)
        for p = 1:length(participant)
            for i = 1:size(GABA_data, 5)
                % devide data at each time point by GMFP
                GABA_data_norm(m, t, p, :, i) = squeeze(GABA_data(m, t, p, 1:30, i)) / GABA_GFP(m, t, p, i);
            end
        end
    end
end
clear m t p i

% calculate DISS and spatial correlation C
for m = 1:length(medication)
    for p = 1:length(participant)
        for i = 1:size(GABA_data_norm, 5)
            diff = squeeze(GABA_data_norm(m, 2, p, :, i) - GABA_data_norm(m, 1, p, :, i));
            GABA_DISS.medication(1, m, p, i) = sqrt(mean(diff.^2));
            GABA_DISS.medication(2, m, p, i) = 1 - (GABA_DISS.medication(1, m, p, i)^2)/2;
        end
    end
end
clear m p i diff

% plot 
for m = 1:length(medication)
    % plot DISS
    data_visual = squeeze(mean(GABA_DISS.medication(:, m, :, :), 3));
    fig = figure(figure_counter); hold on
    plot_DISS(x, data_visual, sprintf('%s: global dissimilarity - %s', target, medication{m}), ...
        colours([(m-1)*2 + 2, (m-1)*2 + 1], :), time_window)
    hold off

    % save the figure
    figure_name = ['DISS_' target '_' medication{m}];
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.png'])

    % update the counter
    figure_counter = figure_counter + 1;
end
clear m data_visual figure_name

% append new variables to the general MATLAB file
save(output_file, 'GABA_data_norm', 'GABA_DISS', '-append');
clear data_temp fig 

%% 7) export data for Ragu
% write text files for Ragu --> uncorrected ppTMS TEPs
for m = 1:length(medication)
    for t = 1:length(time)
        for p = 1:length(participant)
            % choose data to write, remove 'eoi' channels
            data = squeeze(GABA_data(m, t, p, :, :))';

            % define subject
            if participant(p) < 10
                subj = ['S0' num2str(participant(p))];
            else
                subj = ['S' num2str(participant(p))];
            end

            % save as .csv               
            name = ['GABA_' group '_' subj '_' target '_' medication{m} '_' time{t} '.csv']; 
            writematrix(data, [folder_export '\' name])
        end
    end
end
clear m t p data subj name     

% create the montage file
name = [folder_export '\GABA_montage.xyz'];
fileID = fopen(name, 'a');
fprintf(fileID, '30\r\n');
for a = 1:30
    fprintf(fileID, '%.4f %.4f %.4f %s\r\n', ...
        header.chanlocs(a).X, header.chanlocs(a).Y, header.chanlocs(a).Z, header.chanlocs(a).labels);
end
fclose(fileID)
clear name fileID a folder_name

%% 10) TEP amplitude extraction
% ----- decide output parameters -----
GABA_TEP_default.peak = {'P25' 'N40' 'P50' 'N75' 'N100' 'P180'};            % peaks of interest
GABA_TEP_default.center = [0.025 0.04 0.05 0.075 0.10 0.2];                 % default starting latencies
GABA_TEP_default.span = [0.010 0.010 0.010 0.03 0.06 0.06];                 % default peak span   
GABA_TEP_default.eoi = {{'P3' 'CP1' 'CP5'} {'C3' 'FC1' 'FC5'} {'Cz' 'Pz' 'CP1'} {'Fz' 'FC2' 'F4'} {'Fz' 'FC2' 'F4'} {'Cz'}}; 
TEP_polarity = [1 -1 1 -1 -1 1];
percent = 20;                                                               % % of timepoints included in the mean amplitude calculation
map_lims = [-4 4];                                                          % y limits for topoplots -row per stimulus
% ------------------------------------
% set colours for visualisation
col_fig = [1.0000    0.4118    0.1608; 0.9294    0.1412    0.1412; 0.7412    0.0667    0.0667; 
    0.9020    0.1725    0.6588; 0.7176    0.2745    1.0000; 0.3647    0.2078    0.9882];
col_fig1 = [0.0902   0.3725    0.5608; 0.1961    0.5333    0.7608; 0.2549    0.8000    0.8000; 
    0.2549    0.8000    0.5451];

% set n_peaks
n_peaks = 1:length(GABA_TEP_default.peak);   

% loop through subjects and conditions
for p = 1:length(participant)
    % setup names   
    figure_title = sprintf('Subject n. %d', participant(p));    

    % launch summary figure 
    if figure_counter < 3
        figure_counter = 3;
    end
    fig = figure(figure_counter);
    axis_counter = 1;
    
    % adjust figure size
    set(gcf,'units','normalized','outerposition',[0 0 1 1])

    % loop through peaks
    for k = n_peaks
        % identify peak polarity
        polarity = TEP_polarity(k);

        % identify EOIs
        eoi = [];
        for e = 1:length(GABA_TEP_default.eoi{k})
            eoi(e) = find(strcmp(labels, GABA_TEP_default.eoi{k}{e}));
        end

        % choose data
        data_visual = [];
        for m = 1:length(medication)
            for t = 1:length(time)
                % data for timeseries visualization
                data_visual((m-1)*length(medication) + t, :, :) = squeeze(GABA_data(m, t, p, eoi, :)); 
                % data for topoplot
                for e = 1:30
                    data_topoplot((m-1)*length(medication) + t, e, 1, 1, 1, :) = squeeze(GABA_data(m, t, p, e, :));
                end
            end
        end

        % average across eois
        if k == 6
        else
            data_visual = squeeze(mean(data_visual, 2));
        end

        % define default TOI 
        center = GABA_TEP_default.center(k);
        span = GABA_TEP_default.span(k);

        % set manually the final TOI
        finish = 0;
        while finish == 0;                
            % launch the figure
            fig_1 = figure(1);                    

            % plot the background 
            subplot(4, 6, 1:18);
            hold on
            plot(x, data_visual, 'b:', 'LineWidth', 0.5)
            yl = ylim; yl = [yl(1) - 0.2, yl(2) + 0.2]; ylim(yl)
            xlim(time_window)
            rectangle('Position', [-0.005, yl(1), 0.015, yl(2)-yl(1)], 'FaceColor', [0.9020    0.9020    0.9020], 'EdgeColor', 'none')
            title(sprintf('%s\n%s', figure_title, GABA_TEP_default.peak{k}), 'FontSize', 16)
            set(gcf,'units','normalized','outerposition',[0 0 1 1])

            % visualize default peak TOI
            subplot(4, 6, 1:18)
            hold on
            rectangle('Position', [center - span/2, yl(1), span, yl(2)-yl(1)], 'FaceColor', [1,0.7608,0.7608], 'EdgeColor', 'none')

            % calculate mean amplitude
            [y_mean, y_max, lat_peak] = TEP_amplitude(data_visual, center, span, percent, xstep, xstart, polarity);

            % update the figure
            subplot(4, 6, 1:18)
            hold on
            for a = 1:size(data_visual, 1)
                line([lat_peak(a), lat_peak(a)], [0, y_max(a)], 'Color', 'red', 'LineWidth', 1.5)
                plot(x, data_visual(a, :), 'Color', col_fig1(a, :), 'LineWidth', 2.5)
                plot(lat_peak(a), y_max(a), 'o', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'none', 'MarkerSize', 7)
                hold on
            end

            % add lines and params                    
            line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)
            line(x, zeros(length(x)), 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 1)
            set(gca, 'FontSize', 14)
            xlabel('time (s)'); ylabel('GFP (\muV)')

            % add the topoplot   
            for a = 1:size(data_visual, 1)
                subplot(4, 6, 18 + a);
                topo_plot(header, data_topoplot(a, :, :, :, :, :), lat_peak(a), time_window(1), map_lims) 
            end

            % ask for approval
            answer = questdlg('Do you want to proceed?', GABA_TEP_default.peak{k},...
                'Yes, extract outcome values.', 'No, I will adjust TOI manually.', 'Yes, extract outcome values.');

            % switch action
            switch answer
                case 'Yes, extract outcome values.'
                    % close the figure
                    close(fig_1)

                    % exit the while loop
                    finish = 1;

                case 'No, I will adjust TOI manually.'
                    % assign previous center and span
                    choose_center = center;  
                    choose_span = 2 * span;  

                    % identify the limits for visualisation of current peak
                    choose_x1 = ceil((choose_center - choose_span/2 - time_window(1)) / xstep);
                    choose_x2 = ceil((choose_center + choose_span/2 - time_window(1)) / xstep);
                    choose_x = (choose_center - choose_span/2) : xstep : (choose_center + choose_span/2);

                    % prepare data and header for visualization
                    choose_data = data_visual(:, choose_x1 : choose_x2);
                    choose_header = header;
                    choose_header.datasize(6) = length(choose_data);  
                    choose_header.xstart = choose_center - choose_span/2;

                    % check if vector size matches
                    if size(choose_data, 2) ~= length(choose_x)
                        diff = size(choose_data, 2) - length(choose_x);
                        if diff > 0
                            choose_data = choose_data(:, 1:end - diff);
                        elseif diff < 0
                            choose_x = choose_x(1:end + diff);
                        end
                    end

                    % launch the choosing figure                 
                    choose_figure_name = ['Choose manually peak ' GABA_TEP_default.peak{k}];
                    choose_axesHandles = [];
                    choose_fig = figure(2);   
                    choose_axesHandles = [choose_axesHandles subplot(4, 4, [5:16])];  
                    plot(choose_x, choose_data, 'LineWidth', 2)
                    xlim([(choose_center - choose_span/2), (choose_center + choose_span/2)])
                    title(choose_figure_name, 'FontSize', 16)
                    hold on                

                    % plot the line at the center
                    l = line([choose_center, choose_center], get(gca,'ylim'), 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--'); 
                    hold on    

                    % plot the central topography 
                    for a = 1:size(data_visual, 1)
                        choose_axesHandles = [choose_axesHandles subplot(4, 4, a)];
                        topo_plot(header, data_topoplot(a, :, :, :, :, :), choose_center, time_window(1), map_lims);
                    end           

                    % choose the peak position
                    pos_x = get_position(choose_axesHandles);  

                    % update the figure
                    set (choose_fig, 'WindowButtonMotionFcn', '');
                    subplot(4, 4, [5:16])
                    set(l, 'XData', [pos_x, pos_x], 'LineStyle', '-');
                    for a = 1:size(data_visual, 1)
                        subplot(4, 4, a) 
                        cla(choose_axesHandles(2))
                        topo_plot(header, data_topoplot(a, :, :, :, :, :), pos_x, time_window(1), map_lims);
                    end
                    hold off

                    % update the central latency
                    center = pos_x;

                    % close the choosing figure
                    pause(2)
                    close(choose_fig)

                    % close the the main figure
                    close(fig_1)
            end
        end
        clear fig_1 choose_header choose_map_lims choose_span choose_x choose_x1 choose_x2 l a choose_fig pos_x diff...
            choose data choose_center choose_axesHandles answer

        % record outcome variables
        for m = 1:length(medication)
            for t = 1:length(time)
                GABA_TEP_peaks.latency(m, t, p, k) = lat_peak((m-1)*2 + t); 
                GABA_TEP_peaks.amplitude_peak(m, t, p, k) = y_max((m-1)*2 + t); 
                GABA_TEP_peaks.amplitude_mean(m, t, p, k) = y_mean((m-1)*2 + t); 
            end
        end

        % set up the main figure
        figure(fig)
        subplot(length(n_peaks), 7, [1 2 3] + 7*(axis_counter-1))
        hold on
        plot(x, data_visual, ':', 'Color', [0 0.4471 0.7412], 'LineWidth', 0.3)
        yl = get(gca, 'ylim'); 
        xlim(time_window);
        rectangle('Position', [-0.005, yl(1)+0.01, 0.015, yl(2) - yl(1) - 0.02], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')                
        line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 1.5)
        line(x, zeros(length(x)), 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 0.75)
        text(-0.14, 0, sprintf('%s', GABA_TEP_default.peak{k}), 'Color', col_fig(k, :), 'FontSize', 16, 'FontWeight', 'bold')
        set(gca, 'Fontsize', 10)
        ylabel('amplitude (\muV)')
        xlabel('time (s)') 

        % mark peak latencies 
        for a = 1:size(data_visual, 1)
            line([lat_peak(a), lat_peak(a)], [0, y_max(a)], 'Color', col_fig(k, :), 'LineWidth', 2)
            plot(x, data_visual(a, :), 'Color', [0.6 0.6 0.6], 'LineWidth', 0.75)
            hold on
        end 

        % add topoplots
        topoplot_titles = {'placebo - pre' 'placebo - post' 'alprazolam - pre' 'alprazolam - post'};
        for a = 1:size(data_visual, 1)
            subplot(length(n_peaks), 7, 3 + a + 7*(axis_counter-1))
            topo_plot(header, data_topoplot(a, :, :, :, :, :), lat_peak(a), time_window(1), map_lims);   
            if axis_counter == 1
                 set(get(gca, 'title'), 'string', topoplot_titles{a}, 'FontSize', 14);
            end
        end

        % update axis counter
        axis_counter = axis_counter + 1;
    end                                       

    % finalize the summary figure
    figure(fig)  
    subplot(length(n_peaks), 7, [1 2 3])
    hold on     
    yl_sgtitle = get(gca, 'ylim');
    text(-0.14, yl_sgtitle(2)* 1.5, figure_title, 'FontSize', 16, 'FontWeight', 'bold')
    hold off

    % name and save figure
    if participant(p) < 10
        subj = ['0' num2str(participant(p))];
    else
        subj = num2str(participant(p));
    end
    figure_name = ['TEP_' target '_' subj '_amplitude'];
    savefig([folder_figures '\TEP amplitude\' figure_name '.fig'])
    saveas(fig, [folder_figures '\TEP amplitude\' figure_name '.png'])
    close(fig)

    % update the figure counter
    figure_counter = figure_counter + 1;   

% append progressively the output variables to the general MATLAB file
save(output_file, 'GABA_TEP_peaks', '-append');
end
clear p s figure_title k c m t e data_visual data_topoplot fig fig_1 yl center span y_mean y_max lat_peak a ...
    col_fig1 col_fig pos finish axis_counter topoplot_titles eoi n_peaks 

% save data in a R-compatible table 
if ~exist('GABA_TEP_peaks_table_AG')
    GABA_TEP_peaks_table_AG = table;
end
row_counter = height(GABA_TEP_peaks_table_AG) + 1;
for p = 1:length(participant) 
    for m = 1:length(medication)  
        for t = 1:length(time)
            for k = 1:length(GABA_TEP_default.peak) 
                %fill in the table
                GABA_TEP_peaks_table_AG.subject(row_counter) = participant(p);
                GABA_TEP_peaks_table_AG.medication(row_counter) = medication(m);
                GABA_TEP_peaks_table_AG.time(row_counter) = time(t);
                GABA_TEP_peaks_table_AG.stimulus(row_counter) = {'threshold'};
                GABA_TEP_peaks_table_AG.peak(row_counter) = GABA_TEP_default.peak(k);
                GABA_TEP_peaks_table_AG.amplitude_peak(row_counter) = GABA_TEP_peaks.amplitude_peak(m, t, p, k);
                GABA_TEP_peaks_table_AG.amplitude_mean(row_counter) = GABA_TEP_peaks.amplitude_mean(m, t, p, k);
                GABA_TEP_peaks_table_AG.latency(row_counter) = GABA_TEP_peaks.latency(m, t, p, k);

                % update the counter
                row_counter = row_counter + 1;
            end
        end
    end
end
clear m t s p k row_counter
writetable(GABA_TEP_peaks_table_AG, 'GABA_TEP_peaks_table_AG.csv')

% append new variables to the general MATLAB file
save(output_file, 'GABA_TEP_default', 'GABA_TEP_peaks', '-append');
clear percent map_lims

%% 11) plot mean change in TEP peaks
% parameters
col = colours([2, 4], :);
measures = {'mean amplitude' 'peak amplitude' 'latency'};
peak_n = 1:length(GABA_TEP_default.peak);   


% calculate individual TEP change for each peak   
for p = 2:length(participant)
    for m = 1:length(medication)
        for k = 1:length(GABA_TEP_default.peak)
            GABA_TEP_peaks.change(p, m, k, 1) = GABA_TEP_peaks.amplitude_mean(m, 2, p, k)...
                - GABA_TEP_peaks.amplitude_mean(m, 1, p, k);
            GABA_TEP_peaks.change(p, m, k, 2) = GABA_TEP_peaks.amplitude_peak(m, 2, p, k)...
                - GABA_TEP_peaks.amplitude_peak(m, 1, p, k);
            GABA_TEP_peaks.change(p, m, k, 3) = GABA_TEP_peaks.latency(m, 2, p, k)...
                - GABA_TEP_peaks.latency(m, 1, p, k);
        end
    end
end
clear p m k
datasize = size(GABA_TEP_peaks.change);
disp(['Datasize: ' num2str(datasize)])

% calculate group mean values
for m = 1:length(medication)
    for k = 1:length(GABA_TEP_default.peak)
        for a = 1:datasize(end)
            avg_amp(m, k, a) = mean(GABA_TEP_peaks.change(:, m, k, a));
            avg_amp_sem(m, k, a) = std(GABA_TEP_peaks.change(:, m, k, a))/sqrt(length(participant));
        end
    end
end
clear m k a
disp(['Datasize: ' num2str(size(avg_amp))])

% plot a barplot for all peaks together
for a = [2 3]      
    % get data
    data_visual = []; sem_visual = [];
    for m = 1:length(medication)
        peak_counter = 1;
        for k = peak_n
            % identify peak polarity, get the data
            if a ~= 3
                polarity = TEP_polarity(k);
                data_visual(peak_counter, m) = avg_amp(m, k, a) * polarity;
                sem_visual(peak_counter, m) = avg_amp_sem(m, k, a);
            else
                data_visual(peak_counter, m) = avg_amp(m, k, a) * 1000;
                sem_visual(peak_counter, m) = avg_amp_sem(m, k, a) * 1000;
            end

            % update counter
            peak_counter = peak_counter + 1;
        end
    end

    % launch the figure
    fig = figure(figure_counter);
    hold on
    barplot = bar(data_visual, 'EdgeColor', 'none');
    for b = 1:size(data_visual, 2)
        barplot(b).FaceColor = col(b, :);
    end

    % plot errorbars
    ngroups = size(data_visual, 1);
    nbars = size(data_visual, 2);
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    for i = 1:nbars
        x_bar = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x_bar, data_visual(:, i), sem_visual(:, i), sem_visual(:, i), 'k', 'linestyle', 'none');
    end        

    % set other parameters
    title(sprintf('%s: change in TEP %s', target, measures{a}))
    if a ~= 3
        ylabel('amplitude (\muV \pmSEM)');
    else
        ylabel('time (ms \pmSEM)');
    end
    xlabel('TEP component')
    set(gca, 'xtick', 1:length(peak_n), 'xticklabel', GABA_TEP_default.peak(peak_n))
    set(gca, 'Fontsize', 14)
    legend(barplot, medication, 'Location', 'southwest', 'fontsize', 14)

    % name and save figure
    if a ~= 3
        figure_name = ['TEP_' target '_' measures{a}(end - 8:end) '_' measures{a}(1:4) '_all'];
    else
        figure_name = ['TEP_' target '_' measures{a} '_all'];
    end
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.png'])

    % update figure counter
    figure_counter = figure_counter + 1;
end
clear a b m k i peak_n fig fig_name ngroups nbars groupwidth fig barplot ...
    data_visual sem_visual polarity peak_counter x_bar

% plot a boxplot for each TEP peak
for a = 1:size(GABA_TEP_peaks.change, datasize(end))
    for s = [1 2]
        % choose number of displayed peaks
        if s == 1
            peak_n = 2:6;
        else
            peak_n = 1:6;
        end
     
        for k = peak_n 
            % identify peak polarity
            if a ~= 3
                if mod(k, 2) == 0
                    polarity = 1;
                else
                    polarity = -1;
                end
            else
                polarity = 1;
            end
            
            % get the data
            for m = 1:length(medication)
                data_visual(:, m) = GABA_TEP_peaks.change(:, m, s, k, a) * polarity;
            end

            % launch the figure  
            fig = figure(figure_counter);    
            hold on
            boxplot(data_visual, 'color', col)

            % add zero line
            xlim([0.75 2.25])
            line([0.75 2.25], [0, 0], 'LineStyle', ':', 'Color', [0 0 0], 'LineWidth', 1.5)

            % add parameters
            title(sprintf('change in TEP %s:\n%s - peak %s', measures{a}, stimulus{s}, GABA_TEP_default.peak{k}))
            if a ~= 3
                ylabel('amplitude (\muV)');
            else
                ylabel('time (s)');
            end
            xlabel('medication')
            set(gca, 'xtick', 1:length(medication), 'xticklabel', medication)
            set(gca, 'Fontsize', 14)

            % plot the lines
            for p = 1:length(participant)
                p_line(p) = plot([1 2], data_visual(p, [1 2]), '-o',...
                    'Color', [0.75, 0.75, 0.75],...
                    'MarkerSize', 10,...
                    'MArkerEdge', 'none');
                hold on
            end

            % plot the markers
            for m = 1:length(medication)
                scat(m) = scatter(repelem(m, size(data_visual, 1)), data_visual(:, m),...
                    75, col(m, :), 'filled');
            end

            % mark outliers
            h_out = flipud(findobj(gcf,'tag','Outliers'));
            for h = 1:length(h_out)
                x_out =  get(h_out(h), 'XData');
                y_out =  get(h_out(h), 'YData');
                for i = 1:length(x_out)
                    if ~(isnan(x_out(i)))
                        index_out(h, i) = find(data_visual(:, h) == y_out(i));
                        text(x_out(i) + 0.1, double(y_out(i)), sprintf('%d', participant(index_out(h, i))))
                    end
                end
            end
            clear h_out h x_out y_out i index_out

            % name and save figure
            if a ~= 3
                figure_name = ['TEP_' target '_' stimulus{s} '_' measures{a}(end - 8:end) '_' measures{a}(1:4) '_' GABA_TEP_default.peak{k}];
            else
                figure_name = ['TEP_' target '_' stimulus{s} '_' measures{a} '_' GABA_TEP_default.peak{k}];
            end
            savefig([folder_figures '\' figure_name '.fig'])
            saveas(fig, [folder_figures '\' figure_name '.png'])

            % update figure counter
            figure_counter = figure_counter + 1;
        end
    end
end
clear a p k m data_visual scat p_line fig figure_name peak_n polarity

% append new variables to the general MATLAB file
save(output_file, 'GABA_TEP_peaks', '-append');
clear col datasize measures avg_amp avg_amp_sem

%% functions
function peak_x = gfp_plot(x, y, time_window, xstep, labeled, varargin)
% check whether to plot labels (default)
if ~isempty(varargin)
    a = find(strcmpi(varargin, 'max_peaks'));
    if ~isempty(a)
        max_peaks = varargin{a + 1};
    end
end

% launch the figure  
plot(x, y)
yl = get(gca, 'ylim');
cla

% plot interpolated part
hold on
xlim(time_window)
rectangle('Position', [0, yl(1), 0.01, yl(2) - yl(1)], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')

% plot data, mark TMS stimulus
plot(x, y, 'Color', [0 0 0], 'LineWidth', 2.5)
line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)

% find peaks 
[pks, locs] = findpeaks(y, 'MinPeakDistance', 5, 'MinPeakProminence', 0.01);
for a = 1:length(locs)
    if time_window(1) + locs(a)*xstep <= 0.015
        idx(a) = false;
    elseif time_window(1) + locs(a)*xstep > 0.220
        idx(a) = false;        
    else
        idx(a) = true;
    end
end
pks = pks(idx); locs = locs(idx);
if length(pks) > max_peaks
    pks = pks(1:max_peaks); 
    locs = locs(1:max_peaks);
end

% calculate peak coordinations
for a = 1:length(locs)
    peak_x(a) = time_window(1) + locs(a)*xstep;
    peak_y(a) = pks(a);
end
peak_y = double(peak_y);

% plot peaks
for a = 1:length(locs)
    plot(peak_x(a), peak_y(a), ...
        'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none');
    line([peak_x(a), peak_x(a)], [yl(1), peak_y(a)], 'LineStyle', ':', 'Color', [1, 0, 0], 'LineWidth', 1.5)
    
    % label with latency (ms)
    if strcmp(labeled, 'on') 
        text(peak_x(a), peak_y(a) + 0.15, sprintf('%1.0fms', peak_x(a)*1000), 'Color', [1 0 0], 'FontSize', 14)
    end
end

% add parameters
set(gca, 'fontsize', 14)
ylim(yl)
xlabel('time (s)')
ylabel('potential (\muV)')
end
function topo_plot(header, data, x_pos, x_start, map_lims)
varargin = {'maplimits' map_lims 'shading' 'interp' 'whitebk' 'on'};

% fetch data to display
x_visual = ceil((x_pos - x_start)/header.xstep);
vector = data(1, :, 1, 1, 1, x_visual);

%fetch chanlocs
chanlocs = header.chanlocs;

%parse data and chanlocs 
i=1;
for chanpos=1:size(chanlocs,2);
    vector2(i)=double(vector(chanpos));
    chanlocs2(i)=chanlocs(chanpos);
    i=i+1;
end;

topoplot(vector2,chanlocs2,varargin{:});
set(gcf,'color',[1 1 1]);
end
function plot_DISS(x, data_visual, fig_title, colours, time_window)
hold on
title(fig_title, 'FontSize', 16, 'FontWeight', 'bold')
set(gca, 'fontsize', 12); xlabel('time (s)')
ax = gca; ax.XColor = [0.5 0.5 0.5];
plot(x, data_visual(1, :), 'Color', colours(1, :), 'LineWidth', 2.5)
fill([x fliplr(x)],[data_visual(2, :) zeros(1, length(data_visual(2, :)))], ...
    colours(2, :) , 'facealpha', 0.5, 'linestyle', 'none');
line(time_window, [0 0], 'Color', [0, 0, 0], 'LineWidth', 0.5)
line(time_window, [1 1], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5)
line([0, 0], [-0.75 1.75], 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5) 
ylim([-0.75 1.75]); 
lgd = legend({'DISS index' 'spatial correlation'}, 'Location', 'southeast');
lgd.FontSize = 12;
end
function [y_mean, y_max, lat_peak] = TEP_amplitude(data_visual, center, span, percent, xstep, xstart, polarity)
% identify TOI data
start = ceil(((center-span/2) - xstart)/xstep) + 1;
stop = ceil(((center+span/2) - xstart)/xstep);
x_TOI = start : stop;
y_TOI = data_visual(:, x_TOI);

% calculate number of points to average
x_include_n = ceil((percent/100) * size(y_TOI, 2));

% identify peak values for all datasets  
if polarity > 0
    [y_max, x_max] = max(y_TOI, [], 2);
else
    [y_max, x_max] = min(y_TOI, [], 2);
end
x_peak = start + x_max - 1; 

% identify timepoints that will be included in the mean amplitude
x_TOI_include = [];
for a = 1:size(data_visual, 1)
    x_TOI_include_i = x_max(a);
    t_left = ceil((x_include_n - 1)/2); tx_left = x_max(a) - t_left : x_max(a) - 1;
    t_right = x_include_n - 1 - t_left; tx_right = x_max(a) + 1 : x_max(a) + t_right;
    
    % first look left, check for limit
    n = length(find(tx_left <= 0));
    if n > 0
        tx_left = tx_left(tx_left > 0);
        for b = 1:n
            tx_right(end + 1) = tx_right(end) + 1;
        end
    end
    
    % look right, check for limit
    n = length(find(tx_right > size(y_TOI, 2)));
    if n > 0
        tx_right = tx_right(tx_right <= size(y_TOI, 2));
        for b = 1:n
            tx_left(end + 1) = tx_left(1) - b;
        end
    end
    
    % append approved datapoints
    x_TOI_include_i = sort([x_TOI_include_i tx_left tx_right]);
    x_TOI_include = [x_TOI_include; x_TOI_include_i];    
end
x_include = x_TOI_include + start - 1;

% calculate the mean value
for a = 1:size(data_visual, 1)
    y_mean(a, 1) = mean(data_visual(a, x_include(a, :)));
end

% calculate peak latency
lat_peak = x_peak * xstep + xstart;
end
function pos_x = get_position(axesHandles)
% wait until the mouse is clicked
w = waitforbuttonpress;

% get the position of the mouse
CP = get(axesHandles(1), 'CurrentPoint');
pos_x = CP(1,1);

end
