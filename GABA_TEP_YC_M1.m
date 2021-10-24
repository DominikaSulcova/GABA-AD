%% GABA-AD: TMS-EVOKED POTENTIALS - PRIMARY MOTOR CORTEX
% Written by Dominika for GABA-AD project (2020-21)
% 
% % Colection of scripts to analyse preprocessed TMS-EEG data: 
%       * all data and figures are saved in a new folder *
% 1) load the data
%       - loads individual averaged data and trims it in a predefined time
%       window
%       - saves the data from all subjects and conditions in a 6D matrix:
%       2 (medication) X 2 (time) X number of stimuli  X number of subjects X
%       number of channels X number of timepoints
% 
% 2) preliminary visualization of TEPs 
%       - plots baseline TEPs from selected electrode, both conditions in
%       one plot --> mean + CI
% 
% 3) calculate global field power from grand average
%       - averages baseline and post-medication data from both sessions
%       - calculates overall GFP and plots it, adds topoplots for peak
%       times, saves the figure
%       - automatically identifies local peaks within [0.01 0.25]s and
%       saves peak latencies in the outcome structure
% 
% 4) GFP - plot difference
% 
% 5) DISS
% 
% 2) export data for RAGU
%   - removes 'target' channel
%   - saves time-series in a .csv table, timepoint x channel (export folder)
%   - creates an .xyz file with the electrode montage
% 

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
        if exist('colours.mat') > 0
            load('colours.mat')
        else
            disp('No colour scheme found in this directory!')    
        end
end
clear a answer

% create output folders
filename = ['GABA_' group '_' target];
folderpath = uigetdir;
output_file = [folderpath '\' filename '.mat'];

folder_export = [folderpath '\GABA_' group '_export'];
if ~exist(folder_export)
    mkdir(folder_export)
end

folder_figures = [folderpath '\GABA_' group '_figures'];
if ~exist(folder_figures)
    mkdir(folder_figures)  
end

% header
load([prefix ' ' group ' 01 ' target ' ' medication{1} ' ' time{1} ' ' stimulus{3} '.lw6'], '-mat')
labels = {header.chanlocs.labels};

% visualization calculated params
figure_counter = 1;
xstep = header.xstep; 
x = [time_window(1):xstep:time_window(2)];
x_start = (time_window(1) - header.xstart)/xstep;
x_end = (time_window(2) - header.xstart)/xstep;

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
    
                % load individual dataset
                load([prefix ' ' group ' ' subj ' ' target ' ' medication{m} ' ' time{t} ' ' stimulus{s} '.mat'])

                % append the data in the data matrix
                GABA_data(m, t, s, p, :, :) = squeeze(data(:, :, :, :, :, x_start:x_end));                                  
            end
        end
    end
end
clear m t s p data subj
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
    load([folderpath '\' filename '.mat'], 'GABA_data')
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
    load([folderpath '\' filename '.mat'], 'GABA_data')
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

%% 5) DISS
% load data if not loaded
if exist('GABA_GFP') ~= 1
    load([folderpath '\' filename '.mat'], 'GABA_GFP_subject')
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

%% 6) export data for Ragu
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

%% 7) SICI
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

% plot mean SICI GFP - separate graph per medication
for m = 1:length(medication)  
    % launch the figure
    fig = figure(figure_counter);
    hold on

    %  choose data
    data_visual = squeeze(GABA_SICI.GFP_mean(m, :, :));
    CI_visual = squeeze(GABA_SICI.GFP_CI(m, :, :));
    
    % set limits of the figure
    yl = [0 5]; ylim(yl)
    xlim(time_window)

    % shade interpolated interval, add zero line
    plot(x, data_visual(1, :), 'b:', 'LineWidth', 0.5);
    rectangle('Position', [-0.005, yl(1), 0.01, yl(2) - yl(1)], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')
    line(time_window, [0, 0], 'LineStyle', ':', 'Color', [0.75, 0.75, 0.75], 'LineWidth', 1.5)

    % choose colour
    col(1, :) = colours((m-1)*2 + 1, :); col(2, :) = colours((m-1)*2 + 2, :);

    % loop through datasets to plot
    for t = 1:length(time)       
        P(t) = plot(x, data_visual(t, :), 'Color', col(t, :), 'LineWidth', 2.5);
        F(t) = fill([x fliplr(x)],[data_visual(t, :) + CI_visual(t, :) fliplr(data_visual(t, :) - CI_visual(t, :))], ...
            col(t, :), 'FaceAlpha', alpha, 'linestyle', 'none');
    end

    % mark TMS stimulus
    line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)
    


    % add other parameters
    title(sprintf('SICI - GFP: %s', medication{m}))
    xlabel('time (s)')
    ylabel('amplitude (\muV)')
    set(gca, 'FontSize', 14)

    % add legend
    lgd = legend(P, {'baseline' 'post medication'}, 'Location', 'southeast');
    lgd.FontSize = 14;
    hold off
    
    % change figure size
    fig.Position = [250 250 600 350];

    % name and save figure
    figure_name = ['TEP_' target '_SICI_GFP_'  medication{m}];
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.png'])

    % update figure counter
    figure_counter = figure_counter + 1 ;
end
clear e m t s e_n fig sp data_visual CI_visual yl lgd figure_name P F col

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

%% 8) GFP amplitude extraction
% ----- decide output parameters -----
GABA_TEP_default.peak = {'N17' 'P30' 'N45' 'P60' 'N100' 'P180'};            % peaks of interest
GABA_TEP_default.center = [0.017 0.03 0.045 0.06 0.10 0.2];               % default starting latencies
GABA_TEP_default.span = [0.015 0.02 0.015 0.02 0.06 0.06];                  % default peak span
percent = 20;                                                               % % of timepoints included in the mean amplitude calculation
map_lims = [-2.5 2.5; -4 4; -4 4];
% ------------------------------------
% loop through subjects and conditions
for p = 3:length(participant)
    for m = 1:length(medication) 
        for t = 1:length(time)
            for s = 1:length(stimulus)
                % setup names   
                figure_title = sprintf('Subject n. %d: %s, %s medication, %s', ...
                    participant(p), medication{m}, time{t}, stimulus{s});                                       

                % choose data 
                data_visual = double(squeeze(GABA_GFP(m, t, s, p, :))); 
                for e = 1:30
                    data_topoplot(1, e, 1, 1, 1, :) = squeeze(GABA_data(m, t, s, p, e, :));
                end

                % launch summary figure 
                if figure_counter < 3
                    figure_counter = 3;
                end
                fig = figure(figure_counter); 
                hold on

                % initiate the main plot  
                subplot(4, 6, 1:18)
                plot(x, data_visual, ':b')
                yl = get(gca, 'ylim'); 
                xlim(time_window);
                rectangle('Position', [-0.005, yl(1)+0.01, 0.015, yl(2) - yl(1) - 0.02], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')                
                line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)

                % cover for ppTMS
                if t == 1 & s ~= 3
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
                            hold on                      

                            % plot the background 
                            subplot(4, 6, [7:24]);
                            hold on
                            plot(x, data_visual, 'b:', 'LineWidth', 0.5)
                            yl = ylim; yl = [yl(1) - 0.2, yl(2) + 0.2]; ylim(yl)
                            xlim(time_window)
                            rectangle('Position', [-0.005, yl(1), 0.015, yl(2)-yl(1)], 'FaceColor', [0.9020    0.9020    0.9020], 'EdgeColor', 'none')
                            title(sprintf('%s\n%s', figure_title, GABA_TEP_default.peak{k}), 'FontSize', 16)
                            set(gcf,'units','normalized','outerposition',[0 0 1 1])

                            % visualize default peak TOI
                            subplot(4, 6, [7:24])
                            hold on
                            rectangle('Position', [center - span/2, yl(1), span, yl(2)-yl(1)], 'FaceColor', [1,0.7608,0.7608], 'EdgeColor', 'none')

                            % calculate mean amplitude
                            [amplitude, averaged_x, averaged_data] = GFP_amplitude(data_visual', center, span, percent, xstep, time_window(1)); 

                            % calculate  real peak latency (round up)
                            central_latency = averaged_x(ceil(length(averaged_x)/2));

                            % update the figure
                            subplot(4, 6, [7:24])
                            for a = 1:length(averaged_x)
                                line([averaged_x(a), averaged_x(a)], [0, averaged_data(a)], 'Color', 'red', 'LineWidth', 1)
                                hold on
                            end

                            % add the topoplot    
                            subplot(4, 6, 1);
                            topo_plot(header, data_topoplot, central_latency, time_window(1), map_lims(s, :)) 

                            % replot the data to make it visible
                            subplot(4, 6, [7:24])
                            hold on
                            plot(x, data_visual, 'b', 'LineWidth', 2.5)
                            line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)
                            line(x, zeros(length(x)), 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 1)

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
                                    choose_data = double(squeeze(GABA_GFP(m, t, s, p, choose_x1 : choose_x2)));
                                    choose_header = header;
                                    choose_header.datasize(6) = length(choose_data);  
                                    choose_header.xstart = choose_center - choose_span/2;

                                    % check if vector size matches
                                    if length(choose_data) ~= length(choose_x)
                                        diff = length(choose_data) - length(choose_x);
                                        if diff > 0
                                            choose_data = choose_data(1:end - diff);
                                        elseif diff < 0
                                            choose_x = choose_x(1:end + diff);
                                        end
                                    end

                                    % launch the choosing figure                 
                                    choose_figure_name = ['Choose manually peak ' GABA_TEP_default.peak{k}];
                                    choose_axesHandles = [];
                                    choose_fig = figure(2);   
                                    choose_axesHandles = [choose_axesHandles subplot(3, 3, [4:9])];  
                                    plot(choose_x, choose_data, 'LineWidth', 2)
                                    xlim([(choose_center - choose_span/2), (choose_center + choose_span/2)])
                                    title(choose_figure_name, 'FontSize', 16)
                                    hold on                

                                    % plot the line at the center
                                    l = line([choose_center, choose_center], get(gca,'ylim'), 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--'); 
                                    hold on    

                                    % plot the central topography 
                                    choose_axesHandles = [choose_axesHandles subplot(3, 3, 2)];
                                    topo_plot(header, data_topoplot, choose_center, time_window(1), map_lims(s, :));
                                    hold on            

                                    % choose the peak position
                                    pos_x = get_position(choose_axesHandles);  

                                    % update the figure
                                    set (choose_fig, 'WindowButtonMotionFcn', '');
                                    subplot(3, 3, [4:9])
                                    set(l, 'XData', [pos_x, pos_x], 'LineStyle', '-');
                                    subplot(3, 3, 2) 
                                    cla(choose_axesHandles(2))
                                    topo_plot(header, data_topoplot, pos_x, time_window(1), map_lims(s, :));
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
                        clear fig_1 choose_header choose_map_lims choose_span choose_x choose_x1 choose_x2 l

                        % record outcome variables
                        GABA_GFP_peaks.latency(m, t, s, p, k) = central_latency;  
                        GABA_GFP_peaks.amplitude(m, t, s, p, k) = amplitude; 

                        % update the main figure
                        figure(fig)
                        subplot(4, 6, 1:18)
                        hold on

                        %plot the shading
                        for a = 1:length(averaged_x)
                            line([averaged_x(a), averaged_x(a)], [0, averaged_data(a)], 'Color', 'red', 'LineWidth', 1)
                            hold on
                        end

                        % add topoplot
                        subplot(4, 6, 18+k)
                        topo_plot(header, data_topoplot, central_latency, time_window(1), map_lims(s, :));  

                        % shift down
                        pos = get(gca, 'Position');
                        pos(2) = pos(2) - 0.05;
                        set(gca, 'Position', pos);

                        % add timing
                        text(-0.3, -0.8, sprintf('%s', GABA_TEP_default.peak{k}), 'Color', [1 0 0], 'FontSize', 14)
                    end                                       
                else
                    % loop through peaks
                    for k = 1:length(GABA_TEP_default.peak)                 
                        % define default TOI 
                        if t == 2
                            center = GABA_GFP_peaks.latency(m, 1, s, p, k);
                        elseif s == 3
                           center = GABA_GFP_peaks.latency(m, t, 2, p, k); 
                        end
                        span = GABA_TEP_default.span(k);
                        
                        % calculate mean amplitude
                        [amplitude, averaged_x, averaged_data] = GFP_amplitude(data_visual', center, span, percent, xstep, time_window(1)); 

                        % calculate  real peak latency (round up)
                        central_latency = averaged_x(ceil(length(averaged_x)/2));
                        
                        % record outcome variables
                        GABA_GFP_peaks.latency(m, t, s, p, k) = central_latency;  
                        GABA_GFP_peaks.amplitude(m, t, s, p, k) = amplitude; 

                        % update the main figure
                        figure(fig)
                        subplot(4, 6, 1:18)
                        hold on

                        %plot the shading
                        for a = 1:length(averaged_x)
                            line([averaged_x(a), averaged_x(a)], [0, averaged_data(a)], 'Color', 'red', 'LineWidth', 1)
                            hold on
                        end

                        % add topoplot
                        subplot(4, 6, 18+k)
                        topo_plot(header, data_topoplot, central_latency, time_window(1), map_lims(s, :));  

                        % shift down
                        pos = get(gca, 'Position');
                        pos(2) = pos(2) - 0.05;
                        set(gca, 'Position', pos);

                        % add timing
                        text(-0.3, -0.8, sprintf('%s', GABA_TEP_default.peak{k}), 'Color', [1 0 0], 'FontSize', 14)
                    end    
                end
                
                % finalize the summary figure
                figure(fig)
                sgtitle(figure_title)

                % replot the data to make it visible
                subplot(4, 6, 1:18)
                hold on
                plot(x, data_visual, 'color', [0 0 0], 'LineWidth', 2)
                line(x, zeros(length(x)), 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 1)                

                % add other parameters
                set(gca, 'Fontsize', 14)
                ylabel('GFP (\muV)')
                xlabel('time (s)')                
                hold off

                % change figure size
                fig.Position = [250 250 650 550];

                % name and save figure
                if participant(p) < 10
                    subj = ['0' num2str(participant(p))];
                else
                    subj = num2str(participant(p));
                end
                figure_name = ['GFP_' target '_peaks_' subj '_' medication{m} '_' time{t} '_' stimulus{s}];
                savefig([folder_figures '\GFP amplitude\' figure_name '.fig'])
                saveas(fig, [folder_figures '\GFP amplitude\' figure_name '.png'])
                close(fig)

                % update the figure counter
                figure_counter = figure_counter + 1;  
            end
        end
    end    
    % play a celebratory sound at the end of each participant
    tune = load('handel.mat');
    sound(tune.y, tune.Fs)
    
    % append progressively the output variables to the general MATLAB file
    save(output_file, 'GABA_GFP_peaks', '-append');
end
clear p m t s figure_title data_visual e data_topoplot center span fig yl amplitude averaged_x averaged_data central_latency...
    finish answer choose_fig choose_axesHandles choose_center choose_data choose_figure_name pos_x tune   

% save data in a R-compatible table 
if ~exist('GABA_GFP_peak_table')
    GABA_GFP_peak_table = table;
end
row_counter = height(GABA_GFP_peak_table) + 1;
for p = 1:length(participant) 
    for m = 1:length(medication)  
        for t = 1:length(time)
            for s = 1:length(stimulus)
                for k = 1:length(GABA_TEP_default.peak) 
                    %fill in the table
                    GABA_GFP_peak_table.subject(row_counter) = participant(p);
                    GABA_GFP_peak_table.medication(row_counter) = medication(m);
                    GABA_GFP_peak_table.time(row_counter) = time(t);
                    GABA_GFP_peak_table.stimulus(row_counter) = stimulus(s);
                    GABA_GFP_peak_table.peak(row_counter) = GABA_TEP_default.peak(k);
                    GABA_GFP_peak_table.amplitude(row_counter) = GABA_GFP_peaks.amplitude(m, t, s, p, k);
                    GABA_GFP_peak_table.latency(row_counter) = GABA_GFP_peaks.latency(m, t, s, p, k);
                    
                    % update the counter
                    row_counter = row_counter + 1;
                end
            end
        end
    end
end
clear m tm s p k row_counter
writetable(GABA_GFP_peak_table, 'GABA_GFP_peak_table.csv')

% append new variables to the general MATLAB file
save(output_file, 'GABA_TEP_default', 'GABA_GFP_peaks', '-append');

%% GFP amplitude - fixed latency


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
function [amplitude, averaged_x, averaged_data] = GFP_amplitude(data, center, span, percent, step, xstart)
% ------------------------------------------------------------------------
% Fnc: Calculates mean amplitude of the most prominent <percent> of the TOI 
%       --> TOI defined by the center value and the span of the time window
% 
% Input: 
%   data - vector of values (double), GFP 
%   center, span - num values that define the TOI
%   percent - how many top percent of datapoints will be included in the average
%   step, xstart - defines properties of time axes 
% 
% Author: Dominika (2020) 
% ------------------------------------------------------------------------

% prepare the interval vector 
interval = false(1, length(data));

% define the boundaries of the TOI - both limits are included in the interval!
start = round(((center-span/2) - xstart)/step);
stop = round(((center+span/2) - xstart)/step);
data_crop = data(start : stop);

% calculate number of points to average
points_number = ceil((percent/100) * length(data_crop));

% sort the data 
data_sorted = sort(data_crop, 'descend');        

% calculate the mean value
points_included = data_sorted(1:points_number);
amplitude = mean(points_included);

% index averaged points
for i = 1:points_number
    point_pos = find(data_crop == points_included(i)); 
    % in case there are two same values, take the earlier one
    index(i) = point_pos(1);
end
index = sort(index);

% fill ones in the interval vector at indexed positions
for i = 1:points_number
    interval(start + index(i) - 1) = true;
end

% calculate the outcome vectors for visualisation
averaged_data = data .* interval;                                  % keep only values ov averaged datapoints, the rest is set to 0
averaged_index = find(averaged_data);                               % index the averaged datapoints
averaged_x = averaged_index * step + xstart;                        % set the time interval
averaged_data = averaged_data(find(averaged_data));                 % get rid of the zeros

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
