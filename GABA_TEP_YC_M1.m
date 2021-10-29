%% GABA-AD: TMS-EVOKED POTENTIALS - PRIMARY MOTOR CORTEX
% Written by Dominika for GABA-AD project (2020-21)
% 
% Colection of scripts to extract outcome variables from preprocessed TMS-EEG data: 
%   --> extraxted measures are saved in the global output MATLAB file 'GABA_YC_M1_TEPs'
%   --> figures are saved in a folder 'GABA_YC_figures'
%   --> exported datasets (for letswave or Ragu) are saved in a folder 'GABA_YC_export'
% 
% 1) load the data
%       - loads individual averaged data and trims it in a predefined time
%       window
%       - saves the data from all subjects and conditions in a 6D matrix:
%       medication X time X types of stimulus X subjects X channels X timepoints
%       --> saves variable 'GABA_data' to the global output file 
% 
% 2) preliminary visualization of TEPs 
%       - plots baseline TEPs from selected electrode, both conditions in
%       one plot --> mean + CI
%       - plots baseline vs. post-medication 
% 
% 3) calculate individual global field power (GFP)
%       --> saves variable 'GABA_GFP' to the global output file
% 
% 4) identify TEP peaks in grand average GFP
%       - averages baseline and post-medication data from both sessions
%       - calculates overall GFP and plots it, adds topoplots for peak
%       times, saves the figure
%       - automatically identifies local peaks within [0.01 0.25]s and
%       --> saves variable 'GABA_GFP_mean' to the global output file 
% 
% 5) visualize GFP difference
%       - plot baseline vs. post-medication GFP
% 
% 6) DISS
% 
% 7) export data for RAGU
%   - removes 'target' channel
%   - saves time-series in a .csv table, timepoint x channel (export folder)
%   - creates an .xyz file with the electrode montage
% 
% 8) extract GFP amplitude 
% 
% 9) plot mean change in GFP amplitudes 
% 
% 10) extract TEP amplitude 
% 
% 11) plot mean change in TEP peaks
% 
% 12) calculate SICI
% 
% 13) visualize TEP peak amplitude change in SICI
% 
% 14) extract TEP amplitude 
% 
% 15) plot SICI peaks

%% parameters
clear all; clc

% ----- adjustable parameters -----
% dataset
prefix = 'GABA';
group = 'YC';
target = 'M1';
participant = 1:20;
medication = {'placebo' 'alprazolam'}; 
time = {'pre' 'post'};
stimulus = {'CS' 'TS' 'ppTMS'};

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
% load data --> uncorrected ppTMS TEPs
for m = 1:length(medication)
    for t = 1:length(time)
        for s = 1:length(stimulus)
            for p = 1:length(participant)
                % define participant
                if p < 10
                    subj = ['0' num2str(participant(p))];
                else
                    subj = num2str(participant(p));
                end
                
                % identify dataset
                dataset_name = [prefix ' ' group ' ' subj ' ' target ' ' medication{m} ' ' time{t} ' ' stimulus{s}];
                
                % identify cropping limits
                if m == 1 & t == 1 & s == 1 & p == 1
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
                GABA_data(m, t, s, p, :, :) = squeeze(data(:, :, :, :, :, x_start:x_end));                                  
            end
        end
    end
end
clear m t s p data subj dataset_name
disp(['Data import finished. Datasize: ' num2str(size(GABA_data))])

% save dataset to the global MATLAB file
if ~exist() 
    save(output_file, 'GABA_data');
else
    save(output_file, 'GABA_data', '-append');
end

%% 2) preliminary visualization 
% ----- decide output parameters -----
electrode = {'target'};
% ------------------------------------
% load data if not loaded
if exist('GABA_data') ~= 1
    load([folder_results '\' filename '.mat'], 'GABA_data')
end

% prepare average statistics
for m = 1:length(medication)
    for t = 1:length(time)
        for s = 1:length(stimulus)
            for e = 1:size(GABA_data, 5)
                for i = 1:size(GABA_data, 6)
                    GABA_data_mean(m, t, s, e, i) = mean(squeeze(GABA_data(m, t, s, :, e, i)));
                    GABA_data_CI(m, t, s, e, i) = (std(squeeze(GABA_data(m, t, s, :, e, i)))/sqrt(length(participant))) * z;
                end
            end
        end
    end
end
clear m t s e i

% plot the baseline from both sessions
for e = 1:length(electrode) 
    for s = [1 2]
        % identify the electrode
        e_n = find(contains(labels, electrode{e}));
        
        % prepare data
        data_visual = squeeze(GABA_data_mean(:, 1, s, e_n, :));
        CI_visual = squeeze(GABA_data_CI(:, 1, s, e_n, :));

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
        title([target ' - ' stimulus{s} ', ' electrode{e} ' electrode: baseline TEP'])
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
        figure_name = ['TEP_' target '_bl_'  stimulus{s} '_' electrode{e}];
        savefig([folder_figures '\' figure_name '.fig'])
        saveas(fig, [folder_figures '\' figure_name '.png'])

        % update figure counter
        figure_counter = figure_counter + 1 ;
    end
end
clear m s e e_n fig lgd data_visual CI_visual figure_name P F yl

% plot baseline vs. post-medication 
for e = 1:length(electrode) 
    for m = 1:length(medication)
        for s = [1 2]
            % identify the electrode
            e_n = find(contains(labels, electrode{e}));

            % prepare data
            data_visual = squeeze(GABA_data_mean(m, :, s, e_n, :));
            CI_visual = squeeze(GABA_data_CI(m, :, s, e_n, :));

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
            title([target ' - ' stimulus{s} ', ' electrode{e} ' electrode: ' medication{m}])
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
            figure_name = ['TEP_' target '_'  stimulus{s} '_' medication{m} '_' electrode{e}];
            savefig([folder_figures '\' figure_name '.fig'])
            saveas(fig, [folder_figures '\' figure_name '.png'])

            % update figure counter
            figure_counter = figure_counter + 1 ;
        end
    end
end
clear m s t e e_n fig lgd data_visual CI_visual figure_name P F yl

% save dataset to the global MATLAB file
save(output_file, 'GABA_data_mean', 'GABA_data_CI', '-append');
clear electrode 

%% 3) GFP 
% compute individual GFP
for m = 1:length(medication)
    for t = 1:length(time)
        for s = 1:length(stimulus)
            for p = 1:length(participant)
                % calculate GFP (exclude target channel)
                GABA_GFP(m, t, s, p, :) = std(squeeze(GABA_data(m, t, s, p, 1:30, :)), 1);  
            end
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
        for s = 1:length(stimulus)
            GABA_GFP_mean(m, t, s, :) = std(squeeze(GABA_data_mean(m, t, s, 1:30, :)), 1);  
        end
    end 
end
clear m t s

% remove excessive channels in the header
header.chanlocs(31:end) = [];
header.datasize(2) = 30;  

% plot global GFP per stimulus 
for s = 1:length(stimulus)
    % pool all conditions together
    for i = 1:size(GABA_GFP_mean, 4)
        data_visual(i) = mean(GABA_GFP_mean(:, :, s, i), 'all');
    end

    % launch the figure
    fig = figure(figure_counter);
    hold on

    % extract peak latencies
    h_axis(1) = subplot(3, max_peaks, [1 : 2*max_peaks]);
    GABA_peaks = gfp_plot(x, data_visual, time_window, xstep, labeled, 'max_peaks', max_peaks);
    title(['M1 - ' stimulus{s} ', all conditions: GMFP'], 'fontsize', 16, 'fontweight', 'bold')

    % choose data for topoplots 
    for e = 1:30
        for i = 1:size(GABA_data_mean, 5)
            data_topoplot(1, e, 1, 1, 1, i) = mean(GABA_data_mean(:, :, s, e, i), 'all');
        end
    end

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
    hold off

    % save figure
    figure_name = ['GFP_' target '_' stimulus{s}];
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.png'])

    % update figure counteer
    figure_counter = figure_counter + 1;
end
clear s i t e data_visual data_topoplot fig figure_name pos h_axis GABA_peaks

% append new variables to the general MATLAB file
save(output_file, 'GABA_GFP_mean', '-append');
clear labeled max_peaks 

%% 5) GFP - plot difference
% ----- decide output parameters -----
TOI_peaks = [0.017 0.03 0.045 0.060 0.100 0.180];
peaks = {'N17' 'P30' 'N45' 'P60' 'N100' 'P180'};
% ------------------------------------
% plot baseline vs. post-medication
lines = {':' '-'};
for m = 1:length(medication)
    for s = [1 2]
        % load grand average data                
        data_visual = squeeze(mean(squeeze(GABA_GFP(m, :, s, :, :)), 2));

        % launch the figure
        fig = figure(figure_counter);
        h_axis(1) = subplot(4, length(TOI_peaks) + 1 , 1 : 2*(length(TOI_peaks) + 1));
        title(sprintf('GMFP: %s - %s', stimulus{s}, medication{m}),...
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
                for i = 1:size(GABA_data, 6)
                    data_topoplot(1, e, 1, 1, 1, i) = mean(squeeze(GABA_data(m, t, s, :, e, i)), 'all');
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
        figure_name = ['GFP_' target '_' stimulus{s} '_' medication{m}];
        savefig([folder_figures '\' figure_name '.fig'])
        saveas(fig, [folder_figures '\' figure_name '.png'])

        % update figure counteer
        figure_counter = figure_counter + 1;
    end
end
clear peaks lines s m t fig h_axis yl xl F lgd pos k n  P R L descript figure_name e i

%% 6) DISS
% load data if not loaded
if exist('GABA_GFP') ~= 1
    load([folder_results '\' filename '.mat'], 'GABA_GFP_subject')
end

% normalize data by GFP
for m = 1:length(medication)
    for t = 1:length(time)
        for s = 1:length(stimulus) 
            for p = 1:length(participant)
                for i = 1:size(GABA_data, 6)
                    % devide data at each time point by GMFP
                    GABA_data_norm(m, t, s, p, :, i) = squeeze(GABA_data(m, t, s, p, 1:30, i)) / GABA_GFP(m, t, s, p, i);
                end
            end
        end
    end
end
clear m t s p i

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

% ----- baseline -----
% calculate DISS and spatial correlation C
GABA_DISS = struct;
for p = 1:2*length(participant)
    for i = 1:size(GABA_data_norm, 6)
        % between TS and CS stimuli
        diff = squeeze(data_temp(2, p, :, i) - data_temp(1, p, :, i));
        GABA_DISS.baseline.stimulus(1, p, i) = sqrt(mean(diff.^2));
        GABA_DISS.baseline.stimulus(2, p, i) = 1 - (GABA_DISS.baseline.stimulus(1, p, i)^2)/2;
        
        % between ppTMS and TS stimuli
        diff = squeeze(data_temp(3, p, :, i) - data_temp(2, p, :, i));
        GABA_DISS.baseline.SICI(1, p, i) = sqrt(mean(diff.^2));
        GABA_DISS.baseline.SICI(2, p, i) = 1 - (GABA_DISS.baseline.SICI(1, p, i)^2)/2;
    end
end
clear p i diff

% plot baseline TS-CS
data_visual = squeeze(mean(GABA_DISS.baseline.stimulus(:, :, :), 2));
fig = figure(figure_counter); hold on
plot_DISS(x, data_visual, 'Global dissimilarity: baseline TS - CS', colours([2 1], :), time_window)
hold off

% save the figure
figure_name = ['DISS_' target '_bl_stimulus'];
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.png'])

% update the counter
figure_counter = figure_counter + 1;

% plot baseline ppTMS-TS (SICI)
fig = figure(figure_counter); hold on
data_visual = squeeze(mean(GABA_DISS.baseline.SICI(:, :, :), 2));
plot_DISS(x, data_visual, 'Global dissimilarity: baseline SICI', colours([6 5], :), time_window)
hold off

% save the figure
figure_name = ['DISS_' target '_bl_SICI'];
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.png'])

% update the counter
figure_counter = figure_counter + 1;
clear data_visual figure_name

% ----- pre-post medication -----
% calculate DISS and spatial correlation C
for m = 1:length(medication)
    for s = 1:length(stimulus)
        for p = 1:length(participant)
            for i = 1:size(GABA_data_norm, 6)
                diff = squeeze(GABA_data_norm(m, 2, s, p, :, i) - GABA_data_norm(m, 1, s, p, :, i));
                GABA_DISS.medication(1, m, s, p, i) = sqrt(mean(diff.^2));
                GABA_DISS.medication(2, m, s, p, i) = 1 - (GABA_DISS.medication(1, m, s, p, i)^2)/2;
            end
        end
    end
end
clear m s p i diff

% plot 
for m = 1:length(medication)
    for s = 1:length(stimulus) 
        % plot DISS
        data_visual = squeeze(mean(GABA_DISS.medication(:, m, s, :, :), 4));
        fig = figure(figure_counter); hold on
        plot_DISS(x, data_visual, sprintf('Global dissimilarity: %s, %s', medication{m}, stimulus{s}), ...
            colours([(m-1)*2 + 2, (m-1)*2 + 1], :), time_window)
        hold off
        
        % save the figure
        figure_name = ['DISS_' target '_' medication{m} '_' stimulus{s}];
        savefig([folder_figures '\' figure_name '.fig'])
        saveas(fig, [folder_figures '\' figure_name '.png'])
        
        % update the counter
        figure_counter = figure_counter + 1;
    end
end
clear m s data_visual figure_name

% append new variables to the general MATLAB file
save(output_file, 'GABA_data_norm', 'GABA_DISS', '-append');
clear data_temp fig 

%% 7) export data for Ragu
% write text files for Ragu --> uncorrected ppTMS TEPs
for m = 1:length(medication)
    for t = 1:length(time)
        for s = 3
            for p = 1:length(participant)
%                 % choose data to write, remove 'eoi' channels
%                 data = squeeze(GABA_data(m, t, s, p, 1:30, :))';
                
                % define subject
                if participant(p) < 10
                    subj = ['S0' num2str(participant(p))];
                else
                    subj = ['S' num2str(participant(p))];
                end
                
                % load the data
                name_old = ['GABA ' group ' ' subj(2:end) ' ' target ' ' medication{m} ' ' time{t} ' ' stimulus{s} '.mat'];
                load(name_old)
                data = squeeze(data(:, 1:30, :, :, :, x_start:x_end))';
                
                % save as .csv               
                name = ['GABA_' group '_' subj '_' target '_' medication{m} '_' stimulus{s} '_' time{t} '.csv']; 
                writematrix(data, [folder_export '\' name])
                
            end
        end
    end
end
clear m t s p data subj name     

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

%% 8) GFP amplitude extraction
% ----- decide output parameters -----
GABA_TEP_default.peak = {'N17' 'P30' 'N45' 'P60' 'N100' 'P180'};            % peaks of interest
GABA_TEP_default.center = [0.017 0.03 0.045 0.06 0.10 0.2];                 % default starting latencies
GABA_TEP_default.span = [0.010 0.015 0.015 0.03 0.06 0.06];                 % default peak span   
percent = 20;                                                               % % of timepoints included in the mean amplitude calculation
map_lims = [-2.5 2.5; -4 4; -4 4];   
% ------------------------------------
% set colours for visualisation
col_fig = [1.0000    0.4118    0.1608; 0.9294    0.1412    0.1412; 0.7412    0.0667    0.0667; 
    0.9020    0.1725    0.6588; 0.7176    0.2745    1.0000; 0.3647    0.2078    0.9882];
col_fig1 = [0.0902   0.3725    0.5608; 0.1961    0.5333    0.7608; 0.2549    0.8000    0.8000; 
    0.2549    0.8000    0.5451];

% set peak polarity to positive values
polarity = 1;

% loop through subjects and conditions
for p = 4:length(participant)
    for s = 1:length(stimulus)
        % setup names   
        figure_title = sprintf('Subject n. %d: %s stimulus', participant(p), stimulus{s});                                       

        % choose data         
        for m = 1:length(medication)
            for t = 1:length(time)
                % data for timeseries visualization
                data_visual((m-1)*length(medication) + t, :) = squeeze(GABA_GFP(m, t, s, p, :)); 
                % data for topoplot
                for e = 1:30
                    data_topoplot((m-1)*length(medication) + t, e, 1, 1, 1, :) = squeeze(GABA_data(m, t, s, p, e, :));
                end
            end
        end

        % launch summary figure 
        if figure_counter < 3
            figure_counter = 3;
        end
        fig = figure(figure_counter); 
        fig.Position = [2200 500 700 800];

        % initiate the main plot  
        subplot(6, 6, 1:12)
        hold on
        plot(x, data_visual, ':b')
        yl = get(gca, 'ylim'); 
        xlim(time_window);
        rectangle('Position', [-0.005, yl(1)+0.01, 0.015, yl(2) - yl(1) - 0.02], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')                
        line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)

        % loop through peaks
        for k = 1:length(GABA_TEP_default.peak)             
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
                    topo_plot(header, data_topoplot(a, :, :, :, :, :), lat_peak(a), time_window(1), map_lims(s, :)) 
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
                            topo_plot(header, data_topoplot(a, :, :, :, :, :), choose_center, time_window(1), map_lims(s, :));
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
                            topo_plot(header, data_topoplot(a, :, :, :, :, :), pos_x, time_window(1), map_lims(s, :));
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
                choose data choose_center choose_axesHandles answera

            % record outcome variables
            for m = 1:length(medication)
                for t = 1:length(time)
                    GABA_GFP_peaks.latency(m, t, s, p, k) = lat_peak((m-1)*2 + t); 
                    GABA_GFP_peaks.amplitude_peak(m, t, s, p, k) = y_max((m-1)*2 + t); 
                    GABA_GFP_peaks.amplitude_mean(m, t, s, p, k) = y_mean((m-1)*2 + t); 
                end
            end

            % update the main figure
            figure(fig)
            subplot(6, 6, 1:12)
            hold on

            % plot latencies 
            for a = 1:size(data_visual, 1)
                line([lat_peak(a), lat_peak(a)], [0, y_max(a)], 'Color', col_fig(k, :), 'LineWidth', 1.5)
                hold on
            end

            % add topoplots
            for a = 1:size(data_visual, 1)
                subplot(6, 6, 12 + (a-1)*6 + k)
                topo_plot(header, data_topoplot(a, :, :, :, :, :), lat_peak(a), time_window(1), map_lims(s, :)); 
                                
                % shift down
                pos = get(gca, 'Position');
                pos(2) = pos(2) - 0.025;
                set(gca, 'Position', pos);  
                
                % add peak name
                if a == 1                
                    set(get(gca, 'title'), 'string', sprintf('%s', GABA_TEP_default.peak{k}), 'Color', col_fig(k, :), 'FontSize', 14)
                end
            end
        end                                       

        % finalize the summary figure
        figure(fig)
        sgtitle(figure_title)
        hold on

        % replot the data to make it visible
        subplot(6, 6, 1:12)
        hold on
        plot(x, data_visual, 'color', [0.1 0.1 0.1], 'LineWidth', 0.8)              

        % add other parameters
        set(gca, 'Fontsize', 12)
        ylabel('GFP (\muV)')
        xlabel('time (s)')                
        hold off

        % name and save figure
        if participant(p) < 10
            subj = ['0' num2str(participant(p))];
        else
            subj = num2str(participant(p));
        end
        figure_name = ['GFP_' target '_' subj '_amplitude_' '_' stimulus{s}];
        savefig([folder_figures '\GFP amplitude\' figure_name '.fig'])
        saveas(fig, [folder_figures '\GFP amplitude\' figure_name '.png'])
        close(fig)

        % update the figure counter
        figure_counter = figure_counter + 1;  
    end   
    
    % append progressively the output variables to the general MATLAB file
    save(output_file, 'GABA_GFP_peaks', '-append');
end
clear p s figure_title k c m t e data_visual data_topoplot fig fig_1 yl center span y_mean y_max lat_peak a ...
    col_fig1 col_fig pos finish polarity

% save data in a R-compatible table 
if ~exist('GABA_GFP_peaks_table')
    GABA_GFP_peaks_table = table;
end
row_counter = height(GABA_GFP_peaks_table) + 1;
for p = 1:length(participant) 
    for m = 1:length(medication)  
        for t = 1:length(time)
            for s = 1:length(stimulus)
                for k = 1:length(GABA_TEP_default.peak) 
                    %fill in the table
                    GABA_GFP_peaks_table.subject(row_counter) = participant(p);
                    GABA_GFP_peaks_table.medication(row_counter) = medication(m);
                    GABA_GFP_peaks_table.time(row_counter) = time(t);
                    GABA_GFP_peaks_table.stimulus(row_counter) = stimulus(s);
                    GABA_GFP_peaks_table.peak(row_counter) = GABA_TEP_default.peak(k);
                    GABA_GFP_peaks_table.amplitude_peak(row_counter) = GABA_GFP_peaks.amplitude_peak(m, t, s, p, k);
                    GABA_GFP_peaks_table.amplitude_mean(row_counter) = GABA_GFP_peaks.amplitude_mean(m, t, s, p, k);
                    GABA_GFP_peaks_table.latency(row_counter) = GABA_GFP_peaks.latency(m, t, s, p, k);
                    
                    % update the counter
                    row_counter = row_counter + 1;
                end
            end
        end
    end
end
clear m tm s p k row_counter
writetable(GABA_GFP_peaks_table, 'GABA_GFP_peaks_table.csv')

% append new variables to the general MATLAB file
save(output_file, 'GABA_TEP_default', 'GABA_GFP_peaks', '-append');
clear percent map_lims

%% 9) plot mean change in GFP amplitudes  
% calculate indicidual GFP change for each peak   
for p = 2:length(participant)
    for m = 1:length(medication)
        for s = 1:length(stimulus)
            for k = 1:length(GABA_TEP_default.peak)
                GABA_GFP_peaks.amplitude_mean_change(p, m, s, k) = (GABA_GFP_peaks.amplitude_mean(m, 2, s, p, k)...
                    - GABA_GFP_peaks.amplitude_mean(m, 1, s, p, k))/GABA_GFP_peaks.amplitude_mean(m, 1, s, p, k) * 100;
            end
        end
    end
end

% plot as a boxplot
col = colours([2, 4], :);
for s = 1:length(stimulus)
    for k = 1:length(GABA_TEP_default.peak)
    % get the data
    for m = 1:length(medication)
        data_visual(:, m) = GABA_GFP_peaks.amplitude_mean_change(:, m, s, k);
    end

    % launch the figure  
    fig = figure(figure_counter);    
    hold on
    boxplot(data_visual, 'color', col)
    
    % add zero line
    xlim([0.75 2.25])
    line([0.75 2.25], [0, 0], 'LineStyle', ':', 'Color', [0 0 0], 'LineWidth', 1.5)

    % add parameters
    title(sprintf('CHANGE IN GFP: %s - peak %s', stimulus{s}, GABA_TEP_default.peak{k}))
    ylabel('change in GFP (\muV)')
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
    figure_name = ['TEP_' target '_GFP_amplitude_'  stimulus{s} '_' GABA_TEP_default.peak{k}];
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.png'])

    % update figure counter
    figure_counter = figure_counter + 1;
    end
end
clear p s k m data_visual scat col p_line fig figure_name pos_x

%% 10) TEP amplitude extraction
% ----- decide output parameters -----
GABA_TEP_default.eoi.CS = {{} {'Cz' 'CP1' 'C1'} {'Cz' 'FC2' 'C2'} {'C1' 'C3' 'CP1'} {'Cz' 'Fz' 'FC2'} {'Cz'}}; 
GABA_TEP_default.eoi.TS = {{'C3' 'CP1' 'CP5'} {'Cz' 'CP1' 'C1'} {'Cz' 'FC2' 'C2'} {'C1' 'C3' 'CP1'} {'Cz' 'FC1' 'C1'} {'Cz' 'FC2' 'C2'}}; 
percent = 20;                                                               % % of timepoints included in the mean amplitude calculation
map_lims = [-2.5 2.5; -4 4; -4 4];   
% ------------------------------------
% set colours for visualisation
col_fig = [1.0000    0.4118    0.1608; 0.9294    0.1412    0.1412; 0.7412    0.0667    0.0667; 
    0.9020    0.1725    0.6588; 0.7176    0.2745    1.0000; 0.3647    0.2078    0.9882];
col_fig1 = [0.0902   0.3725    0.5608; 0.1961    0.5333    0.7608; 0.2549    0.8000    0.8000; 
    0.2549    0.8000    0.5451];

% loop through subjects and conditions
for p = 2:length(participant)
    for s = 1:length(stimulus)
        % setup names   
        figure_title = sprintf('Subject n. %d: %s stimulus', participant(p), stimulus{s});    
       
        % launch summary figure 
        if figure_counter < 3
            figure_counter = 3;
        end
        fig = figure(figure_counter);
        axis_counter = 1;
        
        % define number of assessed peaks 
        if s == 1
            n_peaks = 2:6;
        else
            n_peaks = 1:6;
        end
        
        % adjust figure size
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        
        % loop through peaks
        for k = n_peaks   
            % identify peak polarity
            if mod(k, 2) == 0
                polarity = 1;
            else
                polarity = -1;
            end
            
            % identify EOIs
            eoi = [];
            if s == 1
                for e = 1:length(GABA_TEP_default.eoi.CS{k})
                    eoi(e) = find(strcmp(labels, GABA_TEP_default.eoi.CS{k}{e}));
                end
            else
                for e = 1:length(GABA_TEP_default.eoi.TS{k})
                    eoi(e) = find(strcmp(labels, GABA_TEP_default.eoi.TS{k}{e}));
                end
            end
            
            % choose data
            data_visual = [];
            for m = 1:length(medication)
                for t = 1:length(time)
                    % data for timeseries visualization
                    data_visual((m-1)*length(medication) + t, :, :) = squeeze(GABA_data(m, t, s, p, eoi, :)); 
                    % data for topoplot
                    for e = 1:30
                        data_topoplot((m-1)*length(medication) + t, e, 1, 1, 1, :) = squeeze(GABA_data(m, t, s, p, e, :));
                    end
                end
            end
            
            % average across eois
            if s == 1 & k == 6
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
                    topo_plot(header, data_topoplot(a, :, :, :, :, :), lat_peak(a), time_window(1), map_lims(s, :)) 
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
                            topo_plot(header, data_topoplot(a, :, :, :, :, :), choose_center, time_window(1), map_lims(s, :));
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
                            topo_plot(header, data_topoplot(a, :, :, :, :, :), pos_x, time_window(1), map_lims(s, :));
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
                    GABA_TEP_peaks.latency(m, t, s, p, k) = lat_peak((m-1)*2 + t); 
                    GABA_TEP_peaks.amplitude_peak(m, t, s, p, k) = y_max((m-1)*2 + t); 
                    GABA_TEP_peaks.amplitude_mean(m, t, s, p, k) = y_mean((m-1)*2 + t); 
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
                topo_plot(header, data_topoplot(a, :, :, :, :, :), lat_peak(a), time_window(1), map_lims(s, :));   
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
        figure_name = ['TEP_' target '_' subj '_amplitude_' '_' stimulus{s}];
        savefig([folder_figures '\TEP amplitude\' figure_name '.fig'])
        saveas(fig, [folder_figures '\TEP amplitude\' figure_name '.png'])
        close(fig)

        % update the figure counter
        figure_counter = figure_counter + 1;  
    end   
    
    % append progressively the output variables to the general MATLAB file
    save(output_file, 'GABA_TEP_peaks', '-append');
end
clear p s figure_title k c m t e data_visual data_topoplot fig fig_1 yl center span y_mean y_max lat_peak a ...
    col_fig1 col_fig pos finish axis_counter topoplot_titles eoi n_peaks 

% save data in a R-compatible table 
if ~exist('GABA_TEP_peaks_table')
    GABA_TEP_peaks_table = table;
end
row_counter = height(GABA_TEP_peaks_table) + 1;
for p = 1:length(participant) 
    for m = 1:length(medication)  
        for t = 1:length(time)
            for s = 1:length(stimulus)
                for k = 1:length(GABA_TEP_default.peak) 
                    %fill in the table
                    GABA_TEP_peaks_table.subject(row_counter) = participant(p);
                    GABA_TEP_peaks_table.medication(row_counter) = medication(m);
                    GABA_TEP_peaks_table.time(row_counter) = time(t);
                    GABA_TEP_peaks_table.stimulus(row_counter) = stimulus(s);
                    GABA_TEP_peaks_table.peak(row_counter) = GABA_TEP_default.peak(k);
                    GABA_TEP_peaks_table.amplitude_peak(row_counter) = GABA_TEP_peaks.amplitude_peak(m, t, s, p, k);
                    GABA_TEP_peaks_table.amplitude_mean(row_counter) = GABA_TEP_peaks.amplitude_mean(m, t, s, p, k);
                    GABA_TEP_peaks_table.latency(row_counter) = GABA_GFP_peaks.latency(m, t, s, p, k);
                    
                    % update the counter
                    row_counter = row_counter + 1;
                end
            end
        end
    end
end
clear m t s p k row_counter
writetable(GABA_TEP_peaks_table, 'GABA_TEP_peaks_table.csv')

% append new variables to the general MATLAB file
save(output_file, 'GABA_TEP_default', 'GABA_TEP_peaks', '-append');
clear percent map_lims

%% 11) plot mean change in TEP peaks
% parameters
col = colours([2, 4], :);
measures = {'mean amplitude' 'peak amplitude' 'latency'};

% calculate individual TEP change for each peak   
for p = 2:length(participant)
    for m = 1:length(medication)
        for s = 1:length(stimulus)
            for k = 1:length(GABA_TEP_default.peak)
                GABA_TEP_peaks.change(p, m, s, k, 1) = GABA_TEP_peaks.amplitude_mean(m, 2, s, p, k)...
                    - GABA_TEP_peaks.amplitude_mean(m, 1, s, p, k);
                GABA_TEP_peaks.change(p, m, s, k, 2) = GABA_TEP_peaks.amplitude_peak(m, 2, s, p, k)...
                    - GABA_TEP_peaks.amplitude_peak(m, 1, s, p, k);
                GABA_TEP_peaks.change(p, m, s, k, 3) = GABA_TEP_peaks.latency(m, 2, s, p, k)...
                    - GABA_TEP_peaks.latency(m, 1, s, p, k);
            end
        end
    end
end
clear p m s k
datasize = size(GABA_TEP_peaks.change);
disp(['Datasize: ' num2str(datasize)])

% calculate group mean values
for m = 1:length(medication)
    for s = 1:length(stimulus)
        for k = 1:length(GABA_TEP_default.peak)
            for a = 1:size(GABA_TEP_peaks.change, datasize(end))
                avg_amp(m, s, k, a) = mean(GABA_TEP_peaks.change(:, m, s, k, a));
                avg_amp_sem(m, s, k, a) = std(GABA_TEP_peaks.change(:, m, s, k, a))/sqrt(length(participant));
            end
        end
    end
end
clear m s k a
disp(['Datasize: ' num2str(size(avg_amp))])

% plot a barplot for all peaks together
for a = 1:size(GABA_TEP_peaks.change, datasize(end))
   for s = [1 2]
        if s == 1
            peak_n = 2:6;
        else
            peak_n = 1:6;
        end     
        
        % get data
        data_visual = [], sem_visual = [];
        for m = 1:length(medication)
            peak_counter = 1;
            for k = peak_n
                % identify peak polarity, get the data
                if a ~= 3
                    if mod(k, 2) == 0
                        polarity = 1;
                    else
                        polarity = -1;
                    end
                    data_visual(peak_counter, m) = avg_amp(m, s, k, a) * polarity;
                    sem_visual(peak_counter, m) = avg_amp_sem(m, s, k, a);
                else
                    data_visual(peak_counter, m) = avg_amp(m, s, k, a) * 1000;
                    sem_visual(peak_counter, m) = avg_amp_sem(m, s, k, a) * 1000;
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
        title(sprintf('change in TEP %s: %s', measures{a}, stimulus{s}))
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
            figure_name = ['TEP_' target '_' stimulus{s} '_' measures{a}(end - 8:end) '_' measures{a}(1:4) '_all'];
        else
            figure_name = ['TEP_' target '_' stimulus{s} '_' measures{a} '_all'];
        end
        savefig([folder_figures '\' figure_name '.fig'])
        saveas(fig, [folder_figures '\' figure_name '.png'])

        % update figure counter
        figure_counter = figure_counter + 1;
    end
end
clear a b s m k i peak_n fig fig_name ngroups nbars groupwidth fig barplot ...
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
clear a p s k m data_visual scat p_line fig figure_name peak_n polarity

% append new variables to the general MATLAB file
save(output_file, 'GABA_TEP_peaks', '-append');
clear col datasize measures avg_amp avg_amp_sem

%% 12) SICI
% ----- decide output parameters -----
electrode = {'C3' 'FC1'};
linestyle = {'-' ':'}
% ------------------------------------
% calculate individual SICI
GABA_SICI = struct;
for m = 1:length(medication)
    for t = 1:length(time)
        for p = 1:length(participant)
            for e = 1:size(GABA_data, 5)
                % calculate GFP (exclude target channel)
                GABA_SICI.individual(m, t, p, e, :) = squeeze(GABA_data(m, t, 3, p, e, :)) - squeeze(GABA_data(m, t, 2, p, e, :));  
            end
        end
    end
end
clear m t p e

% calculate mean SICI
for m = 1:length(medication)
    for t = 1:length(time)
        for e = 1:size(GABA_data, 5)
            for i = 1:size(GABA_data, 6)
                GABA_SICI.mean(m, t, e, i) = mean(squeeze(GABA_SICI.individual(m, t, :, e, i)));
                GABA_SICI.CI(m, t, e, i) = (std(squeeze(GABA_SICI.individual(m, t, :, e, i)))/sqrt(length(participant))) * z;
            end
        end  
    end 
end
clear m t e i 

% plot mean SICI at electrode(s)of interest - baseline, per medication
for m = 1:length(medication) 
    % prepare data
    data_visual = squeeze(GABA_SICI.mean(m, 1, :, :));
    
    % launch the figure
    fig = figure(figure_counter);
    hold on

    % set limits of the figure
    yl = [-5.3 5.3];

    % shade interpolated interval
    plot(x, data_visual(1, :), 'b:', 'LineWidth', 0.5);
    rectangle('Position', [-0.005, yl(1), 0.015, yl(2) - yl(1)], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')

    % choose colour
    col(1, :) = colours((m-1)*2 + 1, :); col(2, :) = colours((m-1)*2 + 2, :);
    
    % plot all channels
    for e = 1:size(data_visual, 1)
        plot(x, data_visual(e, :), 'Color', col(1, :))
    end
    
    % plot electrode(s) of interest
    for e = 1:length(electrode)
        e_n = find(strcmp(labels, electrode{e}));
        e_p(e) = plot(x, data_visual(e_n, :), 'Color', col(2, :), 'LineWidth', 2.5, 'LineStyle', linestyle{e})
    end

    % mark TMS stimulus and zerol line
    line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)
    line(time_window, [0, 0], 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 1.5)

    % add other parameters
    title(sprintf('baseline SICI: %s', medication{m}))
    ylabel('amplitude (\muV)')
    set(gca, 'FontSize', 14)
    xlim(time_window)
    ylim(yl)

    % add legend
    lgd = legend(e_p, electrode, 'Location', 'southeast');
    lgd.FontSize = 14;
    hold off
 
    % change figure size
    fig.Position = [250 250 600 350];

    % name and save figure
    figure_name = ['TEP_' target '_SICI_bl_'  medication{m}];
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.png'])

    % update figure counter
    figure_counter = figure_counter + 1 ;
end
clear m data_visual fig yl col e e_n e_p lgd figure_name linestyle

% plot mean SICI at electrode(s)of interest - per medication
for e = 1:length(electrode) 
    % identify the electrode
    e_n = find(strcmp(labels, electrode{e}));
    
    % loop through conditions
    for m = 1:length(medication)            
        % launch the figure
        fig = figure(figure_counter);
            
        % ----- plot TS and ppTMS separately -----
        for t = 1:length(time)                        
            subplot(3, 1, t) 
            hold on
            
            % prepare data
            data_visual = squeeze(GABA_data_mean(m, t, [2 3], e_n, :));
            CI_visual = squeeze(GABA_data_CI(m, t, [2 3], e_n, :));
            
            % set limits of the figure
            plot(x, data_visual(1, :) + CI_visual(1, :), 'b:', 'LineWidth', 0.5)
            plot(x, data_visual(1, :) - CI_visual(1, :), 'b:', 'LineWidth', 0.5)
            yl = get(gca, 'ylim'); yl(1) = yl(1) - 0.2; yl(2) = yl(2) + 0.3;
            cla, hold on

            % shade interpolated interval 
            plot(x, data_visual(1, :), 'b:', 'LineWidth', 0.5);
            rectangle('Position', [-0.005, yl(1), 0.01, yl(2) - yl(1)], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')
            
            % choose colour
            col(1, :) = colours((m-1)*2 + t, :); col(2, :) = colours(6, :);

            % loop through datasets to plot
            for s = 1:size(data_visual, 1)        
                P(s) = plot(x, data_visual(s, :), 'Color', col(s, :), 'LineWidth', 2.5);
                F(s) = fill([x fliplr(x)],[data_visual(s, :) + CI_visual(s, :) fliplr(data_visual(s, :) - CI_visual(s, :))], ...
                    col(s, :), 'FaceAlpha', alpha, 'linestyle', 'none');
            end
            
            % mark TMS stimulus
            line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)

            % add other parameters
            title([time{t} ' medication'])
            ylabel('amplitude (\muV)')
            set(gca, 'FontSize', 14)
            xlim(time_window)
            ylim(yl)
            
            % add legend
            lgd = legend(P, stimulus([2 3]), 'Location', 'southeast');
            lgd.FontSize = 14;
            hold off
            
            clear P F col
        end
            
    % ----- plot SICI -----
        subplot(3, 1, 3) 
        hold on
        
        data_visual = squeeze(GABA_SICI.mean(m, :, e_n, :));
        CI_visual = squeeze(GABA_SICI.CI(m, :, e_n, :));

        % set limits of the figure
        plot(x, data_visual(1, :) + CI_visual(1, :), 'b:', 'LineWidth', 0.5)
        plot(x, data_visual(1, :) - CI_visual(1, :), 'b:', 'LineWidth', 0.5)
        yl = get(gca, 'ylim'); yl(1) = yl(1) - 0.2; yl(2) = yl(2) + 0.3;
        cla, hold on

        % shade interpolated interval, add zero line
        plot(x, data_visual(1, :), 'b:', 'LineWidth', 0.5);
        rectangle('Position', [-0.005, yl(1), 0.01, yl(2) - yl(1)], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')
        line(time_window, [0, 0], 'LineStyle', ':', 'Color', [0.75, 0.75, 0.75], 'LineWidth', 1.5)

        % choose colour
        col(1, :) = colours((m-1)*2 + 1, :); col(2, :) = colours((m-1)*2 + 2, :);

        % loop through datasets to plot
        for s = 1:size(data_visual, 1)        
            P(s) = plot(x, data_visual(s, :), 'Color', col(s, :), 'LineWidth', 2.5);
            F(s) = fill([x fliplr(x)],[data_visual(s, :) + CI_visual(s, :) fliplr(data_visual(s, :) - CI_visual(s, :))], ...
                col(s, :), 'FaceAlpha', alpha, 'linestyle', 'none');
        end

        % mark TMS stimulus
        line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)

        % add other parameters
        title('SICI')
        xlabel('time (s)')
        ylabel('amplitude (\muV)')
        set(gca, 'FontSize', 14)
        xlim(time_window)
        ylim(yl)

        % add legend
        lgd = legend(P, time, 'Location', 'southeast');
        lgd.FontSize = 14;
        hold off

        clear P F col
            
    % ----- finalize -----
        % add global title
        suptitle(sprintf('SICI: %s, %s electrode', medication{m}, electrode{e}))
        
        % change figure size
        fig.Position = [250 250 500 700];

        % name and save figure
        figure_name = ['TEP_' target '_SICI_'  medication{m} '_' electrode{e}];
        savefig([folder_figures '\' figure_name '.fig'])
        saveas(fig, [folder_figures '\' figure_name '.png'])

        % update figure counter
        figure_counter = figure_counter + 1 ;
        end
end
clear e m t s e_n fig sp data_visual CI_visual yl lgd figure_name

% calculate GFP of SICI
for m = 1:length(medication)
    for t = 1:length(time)
        for p = 1:length(participant)
            % calculate GFP (exclude target channel)
            GABA_SICI.GFP(m, t, p, :) = std(squeeze(GABA_SICI.individual(m, t, p, 1:30, :)), 1);  
        end
    end
end
clear m t p

% calculate mean SICI GFP
for m = 1:length(medication)
    for t = 1:length(time)
        for i = 1:size(GABA_data, 6)
            GABA_SICI.GFP_mean(m, t, i) = mean(squeeze(GABA_SICI.GFP(m, t, :, i)));
            GABA_SICI.GFP_CI(m, t, i) = (std(squeeze(GABA_SICI.GFP(m, t, :, i)))/sqrt(length(participant))) * z;
        end
    end 
end
clear m t i 

% save SICI datasets for letswave
for m = 1:length(medication)
    for t = 1:length(time)
        % create data
        for p = 1:length(participant)
            data(p, :, 1, 1, 1, :) = squeeze(GABA_SICI.individual(m, t, p, :, :));
        end
        
        % create header
        header.name = ['merged SICI ' medication{m} ' ' time{t}];
        header.datasize = size(data);
        header.xstart = time_window(1);
        
        % save
        save([folder_export '\' header.name '.mat'], 'data')
        save([folder_export '\' header.name '.lw6'], 'header')
    end
end
clear m t p 

% append new variable to the general MATLAB file
save(output_file, 'GABA_SICI', '-append');

%% 13) SICI: TEP peak amplitude change
% calculate individual change for peak amplitude of each TEP peak --> only at the baseline  
for p = 1:length(participant)
    for m = 1:length(medication)
        for k = 1:6
            GABA_TEP_peaks.SICI(p, m, k) = GABA_TEP_peaks.amplitude_peak(m, 1, 3, p, k) - GABA_TEP_peaks.amplitude_peak(m, 1, 2, p, k);
        end
    end
end
clear p m k 
        
% get data
data_visual = []; sem_visual = [];
for m = 1:length(medication)
    peak_counter = 1;
    for k = 1:6
        % identify peak polarity, get the data
        if mod(k, 2) == 0
            polarity = 1;
        else
            polarity = -1;
        end
        data_visual(peak_counter, m) = mean(GABA_TEP_peaks.SICI(:, m, k)) * polarity;
        sem_visual(peak_counter, m) = std(GABA_TEP_peaks.SICI(:, m, k))/sqrt(length(participant));

        % update counter
        peak_counter = peak_counter + 1;
    end
end
        
% launch the figure
fig = figure(figure_counter);
hold on
barplot = bar(data_visual, 'EdgeColor', 'none');
col = colours([1, 3], :);
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
title('SICI: change in TEP peak amplitude')
ylabel('amplitude (\muV \pmSEM)');
xlabel('TEP component')
set(gca, 'xtick', 1:6, 'xticklabel', GABA_TEP_default.peak)
set(gca, 'Fontsize', 14)
legend(barplot, medication, 'Location', 'southeast', 'fontsize', 14)

% name and save figure
figure_name = ['TEP_' target '_SICI_amplitude_baseline'];
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.png'])

% update figure counter
figure_counter = figure_counter + 1;
clear a b s m k i peak_n fig fig_name ngroups nbars groupwidth fig barplot ...
    data_visual sem_visual polarity peak_counter x_bar



%% 14) SICI amplitude extraction
% ----- decide output parameters -----
GABA_SICI.default.peak = {'P20' 'N60' 'P200'};           
GABA_SICI.default.center = [0.02 0.06 0.180];                
GABA_SICI.default.span = [0.02 0.04 0.06];               
GABA_SICI.default.eoi = {'C3'}; 
percent = 20;                                                 
map_lims = [-4 4];   
% ------------------------------------
% set colours for visualisation
col_fig = [0.9294    0.1412    0.1412; 0.9020    0.1725    0.6588; 0.7176    0.2745    1];
col_fig1 = [0.0902   0.3725    0.5608; 0.1961    0.5333    0.7608; 0.2549    0.8000    0.8000; 0.2549    0.8000    0.5451];
           
% identify EOI
eoi = find(strcmp(labels, GABA_SICI.default.eoi{1}));

% loop through subjects and conditions
for p = 17:length(participant)
    % setup names   
    figure_title = sprintf('SICI: subject n. %d', participant(p));    

    % launch summary figure 
    if figure_counter < 3
        figure_counter = 3;
    end
    fig = figure(figure_counter);
    axis_counter = 1;

    % adjust figure size
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    
    % choose data
    data_visual = [];
    for m = 1:length(medication)
        for t = 1:length(time)
            % data for timeseries visualization
            data_visual((m-1)*length(medication) + t, :, :) = squeeze(GABA_SICI.individual(m, t, p, eoi, :)); 
            % data for topoplot
            for e = 1:30
                data_topoplot((m-1)*length(medication) + t, e, 1, 1, 1, :) = squeeze(GABA_SICI.individual(m, t, p, e, :));
            end
        end
    end

    % loop through peaks
    for k = 1:length(GABA_SICI.default.peak) 
        % identify peak polarity
        if mod(k, 2) == 0
            polarity = -1;
        else
            polarity = 1;
        end

        % define default TOI 
        center = GABA_SICI.default.center(k);
        span = GABA_SICI.default.span(k);

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
            title(sprintf('%s: %s', figure_title, GABA_SICI.default.peak{k}), 'FontSize', 16)
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
            answer = questdlg('Do you want to proceed?', GABA_SICI.default.peak{k},...
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
                GABA_SICI_peaks.latency(m, t, p, k) = lat_peak((m-1)*2 + t); 
                GABA_SICI_peaks.amplitude_peak(m, t, p, k) = y_max((m-1)*2 + t); 
                GABA_SICI_peaks.amplitude_mean(m, t, p, k) = y_mean((m-1)*2 + t); 
            end
        end

        % set up the main figure
        figure(fig)
        subplot(length(GABA_SICI.default.peak), 7, [1 2 3] + 7*(axis_counter-1))
        hold on
        plot(x, data_visual, ':', 'Color', [0 0.4471 0.7412], 'LineWidth', 0.3)
        yl = get(gca, 'ylim'); 
        xlim(time_window);
        rectangle('Position', [-0.005, yl(1)+0.01, 0.015, yl(2) - yl(1) - 0.02], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')                
        line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 1.5)
        line(x, zeros(length(x)), 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 0.75)
        text(-0.14, 0, sprintf('%s', GABA_SICI.default.peak{k}), 'Color', col_fig(k, :), 'FontSize', 16, 'FontWeight', 'bold')
        set(gca, 'Fontsize', 10)
        ylabel('amplitude (\muV)')
        xlabel('time (s)') 

        % mark peak latencies 
        for a = 1:size(data_visual, 1)
            line([lat_peak(a), lat_peak(a)], [0, y_max(a)], 'Color', col_fig(k, :), 'LineWidth', 2)
            plot(x, data_visual(a, :), 'Color', [0 0.4471 0.7412], 'LineWidth', 0.75)
            hold on
        end 

        % add topoplots
        topoplot_titles = {'placebo - pre' 'placebo - post' 'alprazolam - pre' 'alprazolam - post'};
        for a = 1:size(data_visual, 1)
            subplot(length(GABA_SICI.default.peak), 7, 3 + a + 7*(axis_counter-1))
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
    subplot(length(GABA_SICI.default.peak), 7, [1 2 3])
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
    figure_name = ['TEP_' target '_' subj '_SICI_amplitude'];
    savefig([folder_figures '\SICI amplitude\' figure_name '.fig'])
    saveas(fig, [folder_figures '\SICI amplitude\' figure_name '.png'])
    close(fig)

    % update the figure counter
    figure_counter = figure_counter + 1;  
    
    % append progressively the output variables to the general MATLAB file
    save(output_file, 'GABA_SICI_peaks', '-append');
end
clear p s figure_title k c m t e data_visual data_topoplot fig fig_1 yl center span y_mean y_max lat_peak a ...
    col_fig1 col_fig pos finish axis_counter topoplot_titles eoi n_peaks 

% save data in a R-compatible table 
if ~exist('GABA_SICI_peaks_table')
    GABA_SICI_peaks_table = table;
end
row_counter = height(GABA_SICI_peaks_table) + 1;
for p = 1:length(participant) 
    for m = 1:length(medication)  
        for t = 1:length(time)
            for s = 1:length(stimulus)
                for k = 1:length(GABA_SICI.default.peak) 
                    %fill in the table
                    GABA_SICI_peaks_table.subject(row_counter) = participant(p);
                    GABA_SICI_peaks_table.medication(row_counter) = medication(m);
                    GABA_SICI_peaks_table.time(row_counter) = time(t);
                    GABA_SICI_peaks_table.stimulus(row_counter) = stimulus(s);
                    GABA_SICI_peaks_table.peak(row_counter) = GABA_TEP_default.peak(k);
                    GABA_SICI_peaks_table.amplitude_peak(row_counter) = GABA_SICI_peaks.amplitude_peak(m, t, p, k);
                    GABA_SICI_peaks_table.amplitude_mean(row_counter) = GABA_SICI_peaks.amplitude_mean(m, t, p, k);
                    GABA_SICI_peaks_table.latency(row_counter) = GABA_SICI_peaks.latency(m, t, p, k);
                    
                    % update the counter
                    row_counter = row_counter + 1;
                end
            end
        end
    end
end
clear m t s p k row_counter
writetable(GABA_SICI_peaks_table, 'GABA_SICI_peaks_table.csv')

% append new variables to the general MATLAB file
save(output_file, 'GABA_SICI', 'GABA_SICI_peaks', '-append');
clear percent map_lims

%% 15) plot SICI peaks
% ----- plot baseline -----
% get data
data_visual = []; sem_visual = [];
for m = 1:length(medication)
    for t = 1:length(time)
        peak_counter = 1;
        for k = [1 2]
            data_visual(k, (m-1)*2 + t) = mean(GABA_SICI_peaks.amplitude_peak(m, t, :, k));
            sem_visual(k, (m-1)*2 + t) = std(GABA_SICI_peaks.amplitude_peak(m, t, :, k))/sqrt(length(participant));
        end
    end
end

% launch the figure
fig = figure(figure_counter);
hold on
barplot = bar(data_visual, 'EdgeColor', 'none');
for b = 1:size(data_visual, 2)
    barplot(b).FaceColor = colours(b, :);
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
title('change in SICI peak amplitude')
ylabel('amplitude (\muV \pmSEM)');
xlabel('SICI component')
set(gca, 'xtick', [1 2], 'xticklabel', GABA_SICI.default.peak([1 2]))
set(gca, 'Fontsize', 14)
legend(barplot, {'placebo - pre' 'placebo - post' 'alprazolam - pre' 'alprazolam - post'}, 'Location', 'northeast', 'fontsize', 14)

% name and save figure
figure_name = ['TEP_' target '_SICI_amplitude_all'];
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.png'])

% update figure counter
figure_counter = figure_counter + 1;
clear b s m k i t peak_n fig fig_name ngroups nbars groupwidth fig barplot ...
    data_visual sem_visual polarity peak_counter x_bar

% ----- change introduced by medication -----
% calculate individual SICI change for each peak   
for p = 1:length(participant)
    for m = 1:length(medication)
        for k = [1 2]
            GABA_SICI_peaks.change(p, m, k, 1) = GABA_SICI_peaks.amplitude_mean(m, 2, p, k)...
                - GABA_SICI_peaks.amplitude_mean(m, 1, p, k);
            GABA_SICI_peaks.change(p, m, k, 2) = GABA_SICI_peaks.amplitude_peak(m, 2, p, k)...
                - GABA_SICI_peaks.amplitude_peak(m, 1, p, k);
            GABA_SICI_peaks.change(p, m, k, 3) = GABA_SICI_peaks.latency(m, 2, p, k)...
                - GABA_SICI_peaks.latency(m, 1, p, k);
        end
    end
end
clear p m s k
datasize = size(GABA_SICI_peaks.change);
disp(['Datasize: ' num2str(datasize)])

% calculate group mean change of peak amplitude
for m = 1:length(medication)
    for k = [1 2]
        avg_amp(m, k) = mean(GABA_SICI_peaks.change(:, m, k, 2));
        avg_amp_sem(m, k) = std(GABA_TEP_peaks.change(:, m, k, 2))/sqrt(length(participant));
    end
end
clear m s k a
disp(['Datasize: ' num2str(size(avg_amp))])
 
% plot a boxplot for each SICI peak   
col = colours([2, 4], :);
for k = [1 2]
    % identify peak polarity
    if mod(k, 2) == 0
        polarity = -1;
    else
        polarity = 1;
    end

    % get the data
    for m = 1:length(medication)
        data_visual(:, m) = GABA_SICI_peaks.change(:, m, k, 2) * polarity;
    end

    % launch the figure  
    fig = figure(figure_counter);    
    hold on
    boxplot(data_visual, 'color', col)

    % add zero line
    xlim([0.75 2.25])
    line([0.75 2.25], [0, 0], 'LineStyle', ':', 'Color', [0 0 0], 'LineWidth', 1.5)

    % add parameters
    title(sprintf('change in SICI peak amplitude: %s', GABA_SICI.default.peak{k}))
    ylabel('amplitude (\muV)');
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
    figure_name = ['TEP_' target '_SICI_amplitude_' GABA_SICI.default.peak{k}];
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.png'])

    % update figure counter
    figure_counter = figure_counter + 1;
end
clear a p k m data_visual scat p_line fig figure_name peak_n polarity

% append new variables to the general MATLAB file
save(output_file, 'GABA_SICI_peaks', '-append');
clear col subj datasize measures avg_amp avg_amp_sem

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
function mouse_move(hObject,eventdata, axesHandles, header_visual, data_visual)

% get the position of the mouse
CP = get(hObject, 'CurrentPoint');
position = get(hObject, 'Position');
xCursor = CP(1,1)/position(1,3); % normalize
yCursor = CP(1,2)/position(1,4); % normalize

% get the position of the axes within the GUI
axesPos = get(axesHandles(1),'Position');
minx    = axesPos(1);
miny    = axesPos(2);
maxx    = minx + axesPos(3);
maxy    = miny + axesPos(4);

% check if the mouse is within the axes 
if xCursor >= minx && xCursor <= maxx && yCursor >= miny && yCursor <= maxy
    % get the cursor position within the lower axes
    currentPoint = get(axesHandles(1),'CurrentPoint');
    x = currentPoint(2,1);
    y = currentPoint(2,2);
    % update the topoplot in the uper axes
    map_lims = get(axesHandles(1), 'ylim');
    cla(axesHandles(2))
    subplot(3, 3, 2)
    TEP_topoplot(header_visual, data_visual, x, map_lims); 
    hold on
end
end
function pos_x = get_position(axesHandles)
% wait until the mouse is clicked
w = waitforbuttonpress;

% get the position of the mouse
CP = get(axesHandles(1), 'CurrentPoint');
pos_x = CP(1,1);

end
function TEP_topoplot(header, data, x_pos, map_lims)
varargin = {'maplimits' map_lims 'shading' 'interp' 'whitebk' 'on'};

% fetch data to display
x_visual = ceil((x_pos - header.xstart)/header.xstep);
vector = squeeze(data(1,:,1,1,1,x_visual));

%fetch chanlocs
chanlocs = header.chanlocs;

%parse data and chanlocs 
i = 1;
for chanpos=1:size(chanlocs,2);
    vector2(i)=double(vector(chanpos));
    chanlocs2(i)=chanlocs(chanpos);
    i = i + 1;
end;

topo_plot(vector2,chanlocs2,varargin{:});
set(gcf,'color',[1 1 1]);
end
