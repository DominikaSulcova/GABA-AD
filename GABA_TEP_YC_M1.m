%% GABA-AD: TMS-EVOKED POTENTIALS
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
% 3) calculate global mean field power from grand average
%       - averages baseline and post-medication data from both sessions
%       - calculates overall GMFP and plots it, adds topoplots for peak
%       times, saves the figure
%       - automatically identifies local peaks within [0.01 0.25]s and
%       saves peak latencies in the outcome structure
% 
% 4) define peak widths using grand average GMFP
%       - for each stimulus separately, plots and saves the figure 
%       - uses function findpeaks awith the 'halfheight' option
%       - saves peak widths in the outcome structure
% 
% 5) determine electrodes of interest
%       - uses GMFP peak latencies as default search points and signal from Cz 
%       as starting peak point to track peak progress
%       - averages across subjects
%       - identifies 3 electrodes with the highest average amplitude for
%       each peak, plots the  




%% parameters
clear all; clc

% ----- adjustable parameters -----
% dataset
group = 'YC';
target = 'M1';
participant = 1:20;
medication = {'placebo' 'alprazolam'}; 
time = {'pre' 'post'};

% visualization 
time_window = [-0.05, 0.3];
z = 1.96;
alpha = 0.2;
% --------------------------------

% dataset parameters
stimulus = {'CS' 'TS' 'ppTMS'};
prefix = 'eois crop avg final_dataset';
load([prefix ' ' stimulus{1} ' ' group num2str(participant(1)) ' ' medication{1} ' ' time{1} ' ' target '.lw6'], '-mat')
labels.CS = {header.chanlocs.labels};
load([prefix ' ' stimulus{2} ' ' group num2str(participant(1)) ' ' medication{1} ' ' time{1} ' ' target '.lw6'], '-mat')
labels.TS = {header.chanlocs.labels};

% visualization calculated params
figure_counter = 1;
xstep = header.xstep; 
x = [time_window(1):xstep:time_window(2)];
x_start = (time_window(1) - header.xstart)/xstep;
x_end = (time_window(2) - header.xstart)/xstep;

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
folderpath = [pwd '\' filename];
if ~exist(folderpath) 
    mkdir(folderpath)
    mkdir(folderpath, [filename '_figures'])
end      

%% 1) extract individual data
% load data --> corrected ppTMS TEPs
n = length(labels.TS) - length(labels.CS);
for m = 1:length(medication)
    for t = 1:length(time)
        for s = 1:length(stimulus)
            for p = 1:length(participant)
                % load individual dataset
                load([prefix ' ' stimulus{s} ' ' group num2str(participant(p)) ' ' medication{m} ' ' time{t} ' ' target '.mat'])

                % append the data in the data matrix
                GABA_data(m, t, s, p, [1:size(data, 2)], :) = squeeze(data(:, :, :, :, :, x_start:x_end));   

                % fill in zeros if n of electrodes doesn't match 
                if n > 0 & strcmp(stimulus{s}, 'CS') | n < 0 & strcmp(stimulus{s}, 'TS')
                    GABA_data(m, t, s, p, [size(data, 2)+1: size(data, 2)+n], :) = zeros(n, size(GABA_data, 6));
                end                                
            end
        end
    end
end
clear m t s p data n
disp(['Data import finished. Datasize: ' num2str(size(GABA_data))])

% save dataset to the global MATLAB file
if ~exist([folderpath '\' filename '.mat']) 
    save([folderpath '\' filename '.mat'], 'GABA_data');
else
    save([folderpath '\' filename '.mat'], 'GABA_data', '-append');
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
                    data_mean(m, t, s, e, i) = mean(squeeze(GABA_data(m, t, s, :, e, i)));
                    data_CI(m, t, s, e, i) = (std(squeeze(GABA_data(m, t, s, :, e, i)))/sqrt(length(participant))) * z;
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
        statement = ['e_n = find(contains(labels.' stimulus{s} ', electrode{e}));'];
        eval(statement)
        
        % prepare data
        data_visual = squeeze(data_mean(:, 1, s, e_n, :));
        CI_visual = squeeze(data_CI(:, 1, s, e_n, :));

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
        savefig([folderpath '\' filename '_figures\' figure_name '.fig'])
        saveas(fig, [folderpath '\' filename '_figures\' figure_name '.png'])

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
            statement = ['e_n = find(contains(labels.' stimulus{s} ', electrode{e}));'];
            eval(statement)

            % prepare data
            data_visual = squeeze(data_mean(m, :, s, e_n, :));
            CI_visual = squeeze(data_CI(m, :, s, e_n, :));

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
            savefig([folderpath '\' filename '_figures\' figure_name '.fig'])
            saveas(fig, [folderpath '\' filename '_figures\' figure_name '.png'])

            % update figure counter
            figure_counter = figure_counter + 1 ;
        end
    end
end
clear m s t e e_n fig lgd data_visual CI_visual figure_name P F yl

% save dataset to the global MATLAB file
save([folderpath '\' filename '.mat'], 'data_mean', 'data_CI', '-append');
clear electrode 

%% 3) GMFP - peak identification
% ----- decide output parameters -----
labeled = 'off';
max_peaks = 6;
% ------------------------------------
if exist('GABA_data') ~= 1
    load([folderpath '\' filename '.mat'], 'GABA_data')
end

% calculate GMFP (exclude target and eoi channels)
for m = 1:length(medication)
    for t = 1:length(time)
        for s = 1:length(stimulus)
            GABA_GMFP(m, t, s, :) = std(squeeze(data_mean(m, t, s, 1:30, :)), 1);  
        end
    end 
end
clear m t s

% remove excessive channels in the header
header.chanlocs(31:end) = [];
header.datasize(2) = 30;  

% plot global GMFP per stimulus 
for s = 1:length(stimulus)
    % pool all conditions together
    for i = 1:size(GABA_GMFP, 4)
        data_visual(i) = mean(GABA_GMFP(:, :, s, i), 'all');
    end

    % launch the figure
    fig = figure(figure_counter);
    hold on

    % extract peak latencies
    h_axis(1) = subplot(3, max_peaks, [1 : 2*max_peaks]);
    GABA_peaks = gmfp_plot(x, data_visual, time_window, xstep, labeled, 'max_peaks', max_peaks);
    title(['M1 - ' stimulus{s} ', all conditions: GMFP'], 'fontsize', 16, 'fontweight', 'bold')

    % choose data for topoplots 
    for e = 1:30
        for i = 1:size(data_mean, 5)
            data_topoplot(1, e, 1, 1, 1, i) = mean(data_mean(:, :, s, e, i), 'all');
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
    figure_name = ['GMFP_' target '_' stimulus{s}];
    savefig([folderpath '\' filename '_figures\' figure_name '.fig'])
    saveas(fig, [folderpath '\' filename '_figures\' figure_name '.png'])

    % update figure counteer
    figure_counter = figure_counter + 1;
end
clear s i t e data_visual data_topoplot fig figure_name pos h_axis GABA_peaks

% append new variables to the general MATLAB file
save([folderpath '\' filename '.mat'], 'GABA_GMFP', '-append');
clear labeled max_peaks 

%% 4) GMFP - plot difference
% ----- decide output parameters -----
TOI_peaks = [0.03 0.045 0.060 0.100 0.180];
peaks = {'P30' 'N45' 'P60' 'N100' 'P180'};
% ------------------------------------
% calculate GMFP for each subject/condition
for m = 1:length(medication)
    for t = 1:length(time)
        for s = 1:length(stimulus) 
            for p = 1:length(participant)
                % calculate GMFP
                GABA_GMFP_subject(m, t, s, p, :) = std(squeeze(GABA_data(m, t, s, p, 1:30, :)), 1);
            end
        end
    end
end
clear m t s p 

% plot baseline vs. post-medication
lines = {':' '-'};
for m = 1:length(medication)
    for s = [1 2]
        % load grand average data                
        data_visual = squeeze(mean(squeeze(GABA_GMFP_subject(m, :, s, :, :)), 2));

        % launch the figure
        fig = figure(figure_counter);
        h_axis(1) = subplot(4, length(TOI_peaks) + 1 , 1 : 2*(length(TOI_peaks) + 1));
        title(sprintf('GMFP: %s - %s', stimulus{s}, medication{m}),...
            'FontSize', 16, 'FontWeight', 'bold')
        set(gca, 'fontsize', 12)
        xlabel('time (s)')
        ylabel('GMFP (\muV)')
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

        % plot GMFP
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
        figure_name = ['GMFP_' target '_' stimulus{s} '_' medication{m}];
        savefig([folderpath '\' filename '_figures\' figure_name '.fig'])
        saveas(fig, [folderpath '\' filename '_figures\' figure_name '.png'])

        % update figure counteer
        figure_counter = figure_counter + 1;
    end
end
clear lines s m t fig h_axis yl xl F lgd pos k n  P R L descript figure_name e i

% append new variables to the general MATLAB file
save([folderpath '\' filename '.mat'], 'GABA_GMFP_subject', '-append');
clear TOI TOI_peaks peaks

%% 5) DISS
% load data if not loaded
if exist('GABA_GMFP_subject') ~= 1
    load([folderpath '\' filename '.mat'], 'GABA_GMFP_subject')
end

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
savefig([folderpath '\' filename '_figures\' figure_name '.fig'])
saveas(fig, [folderpath '\' filename '_figures\' figure_name '.png'])

% update the counter
figure_counter = figure_counter + 1;

% plot baseline ppTMS-TS (SICI)
fig = figure(figure_counter); hold on
data_visual = squeeze(mean(GABA_DISS.baseline.SICI(:, :, :), 2));
plot_DISS(x, data_visual, 'Global dissimilarity: baseline SICI', colours([6 5], :), time_window)
hold off

% save the figure
figure_name = ['DISS_' target '_bl_SICI'];
savefig([folderpath '\' filename '_figures\' figure_name '.fig'])
saveas(fig, [folderpath '\' filename '_figures\' figure_name '.png'])

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
        savefig([folderpath '\' filename '_figures\' figure_name '.fig'])
        saveas(fig, [folderpath '\' filename '_figures\' figure_name '.png'])
        
        % update the counter
        figure_counter = figure_counter + 1;
    end
end
clear m s data_visual figure_name

% append new variables to the general MATLAB file
save([folderpath '\' filename '.mat'], 'GABA_data_norm', 'GABA_DISS', '-append');
clear data_temp fig 

%% 6) export data for Ragu
% ----- adjustable parameters -----
folder_name = {'corrected' 'not_corrected'};
% ----- adjustable parameters -----

% create output folder
for a = 1:length(folder_name)
    mkdir([folderpath '\' filename '_export\' folder_name{a}])
end
clear a

% write text files for Ragu --> corrected ppTMS TEPs
for m = 1:length(medication)
    for t = 1:length(time)
        for s = 1:length(stimulus)
            for p = 1:length(participant)
                % choose data to write, remove 'eoi' channels
                data = squeeze(GABA_data(m, t, s, p, 1:30, :))';
                
                % define subject
                if participant(p) < 10
                    subj = ['S0' num2str(participant(p))];
                else
                    subj = ['S' num2str(participant(p))];
                end
                
                % save as .csv               
                name = ['GABA_' subj '_' medication{m} '_' stimulus{s} '_' time{t} '.csv']; 
                writematrix(data, [folderpath '\' filename '_export\' folder_name{1} '\' name])
            end
        end
    end
end
clear m t s p data subj name     

% % write text files for Ragu --> uncorrected ppTMS TEPs
% for m = 1:length(medication)
%     for t = 1:length(time)
%         for s = 1:length(stimulus)
%             for p = 1:length(participant)
%                 % choose data to write, remove 'eoi' channels
%                 data = squeeze(GABA_data(m, t, s, p, 1:30, :))';
%                 
%                 % define subject
%                 if participant(p) < 10
%                     subj = ['S0' num2str(participant(p))];
%                 else
%                     subj = ['S' num2str(participant(p))];
%                 end
%                 
%                 % save as .csv               
%                 name = ['GABA_' subj '_' medication{m} '_' stimulus{s} '_' time{t} '.csv']; 
%                 writematrix(data, [folderpath '\' filename '_export\' folder_name{2} '\' name])
%             end
%         end
%     end
% end
% clear m t s p data subj name   

% create the montage file
name = [folderpath '\' filename '_export\GABA_montage.xyz'];
fileID = fopen(name, 'a');
fprintf(fileID, '30\r\n');
for a = 1:30
    fprintf(fileID, '%.4f %.4f %.4f %s\r\n', ...
        header.chanlocs(a).X, header.chanlocs(a).Y, header.chanlocs(a).Z, header.chanlocs(a).labels);
end
fclose(fileID)
clear name fileID a folder_name

%% ) peak widths
% calculate narrow window parameters
x_start_narrow = (0.01 - time_window(1))/xstep;
x_end_narrow = (0.25 - time_window(1))/xstep;
x_narrow = [0.01:xstep:0.25];

% loop through stimuli
for s = 1:length(stimulus)
    % average data across conditions
    for i = 1:length(x_narrow)
        data_visual(i) = mean(GABA_GMFP(:, :, s, x_start_narrow + (i - 1)), 'all');
    end
        
    % identify peak widths
    [P, L, GABA_peaks(s).widths, R] = findpeaks(data_visual, 'Annotate','extents', 'WidthReference', 'halfheight');
    GABA_peaks(s).widths = ceil(GABA_peaks(s).widths)* header.xstep;
    
    % plot figure
    fig = figure(figure_counter);
    hold on
    findpeaks(data_visual, x_narrow, 'Annotate', 'extents', 'WidthReference', 'halfheight');
    set(gca, 'fontsize', 14)
    xlabel('time(s)')
    ylabel('power (\muV^2)')
    grid off
    
    % add width denotation
    for k = 1:length(GABA_peaks(s).widths)
        if k == length(GABA_peaks(s).widths)
            text(GABA_peaks(s).latencies(k) - 0.005, -0.25, sprintf('%1.0f ms', GABA_peaks(s).widths(k)*1000), 'Color', [0.93 0.69 0.13], 'FontSize', 14)
        else
            text(GABA_peaks(s).latencies(k) - 0.005, -0.25, sprintf('%1.0f', GABA_peaks(s).widths(k)*1000), 'Color', [0.93 0.69 0.13], 'FontSize', 14)
        end
    end
    xlim([-0.005, 0.26])
    
    % save figure, update
    savefig(['GABA_' target '_widths'])
    saveas(fig, ['GABA_' target '_widths.png'])
    figure_counter = figure_counter + 1;
    end
clear s i

%%% ONLY FOR AG %%%
GABA_peaks(s).widths = [GABA_peaks(s).widths(1:2), GABA_peaks(s).widths(2), GABA_peaks(s).widths(3:5)];
%%%%%%%%%%%%%%%%%%%

% append new variables to the general MATLAB file
save(filename, 'GABA_TEP', '-append');
clear x_start_narrow x_end_narrow x_narrow s k data fig P L R 

%% 5) EOIs
% ----- decide output parameters -----
seed_electrode = {'target' 'Cz'};                           % electrode that will be used to set 0 timepoint
seed_peaks = {1:4, 5:6};                                    % which peek s use which seed electrode
buffer = 0.5;                                               % a margin of the window for peak visualisation 
eoi_n = 3;                                                  % number of detected EOIs
% ------------------------------------

% check if numbers match
if numel(seed_electrode) ~= numel(seed_peaks)
    disp('Number of seed electrodes does not correspond to the seed peak distribution!')
end

% k = 1; p = 1; m = 1; s = 1; e = 1;
for k = 1:length(GABA_peaks(1).peaks)
    % select seed 
    for a = 1:numel(seed_electrode)
        if any(seed_peaks{a} == k)
            seed = find(contains(labels, seed_electrode{a}));
        end
    end
    
    % loop through the datasets
    counter = 1;
    for p = 1:length(participant)
        for m = 1:length(medication)
            for s = 1:length(stimulus)
                % choose data from the seed electrode
                for e = 1:size(GABA_data, 5)
                    data(1, e, 1, 1, 1, :) = squeeze(GABA_data(m, 1, s, p, e, :));
                end

                % identify peak latency for current dataset
                [peak_center, data_c, data_c_sub] = track_peak(data, header, time_window, k, GABA_peaks(s), buffer, seed);
                
                % fill in outcome structure
                GABA_TEP_subject(p).latency(m, k) = peak_center;
                
                % append centered data
                statement = ['GABA_' GABA_peaks(s).peaks{k} '_centered(counter, :, :) = data_c;'];
                eval(statement)
                
                % append centered subtracted data
                statement = ['GABA_' GABA_peaks(s).peaks{k} '_subtracted(counter, :, :) = data_c_sub;'];
                eval(statement)
                
                % update row counter
                counter = counter + 1;
            end
        end
    end
end
clear counter a k p m t s e data statement peak_center data_c data_c_sub seed

% average across subjects and sessions, save for letswave
for k = 1:length(GABA_peaks(1).peaks)
    % choose centered data
    statement = ['data_i =  GABA_' GABA_peaks(1).peaks{k} '_centered;'];
    eval(statement)
    
    % average across trials (subject x medication)
    for e = 1:size(data_i, 2)
        for i = 1:size(data_i, 3)
            data(1, e, 1, 1, 1, i) =  mean(squeeze(data_i(:, e, i)));         
            GABA_tracking(k).centered(e, i) = mean(squeeze(data_i(:, e, i)));        
        end
    end
    
    % save for LW
    filename = ['TEP tracked ' target ' '  GABA_peaks(1).peaks{k} ' centered'];
    header.name = filename; 
    header.datasize(6) = size(data, 6);
    span = (1 + buffer) * GABA_peaks(1).widths(k);
    header.xstart = - span/2;
    save([filename '.mat'], 'data');
    save([filename '.lw6'], 'header');    
    clear data
    
    % choose centered subtracted data
    statement = ['data_i =  GABA_' GABA_peaks(1).peaks{k} '_subtracted;'];
    eval(statement)
    
    % average across trials (subject x medication)
    for e = 1:size(data_i, 2)
        for i = 1:size(data_i, 3)
            data(1, e, 1, 1, 1, i) =  mean(squeeze(data_i(:, e, i)));         
            GABA_tracking(k).subtracted(e, i) = mean(squeeze(data_i(:, e, i)));        
        end
    end
    
    % save for LW
    filename = ['TEP tracked ' target ' '  GABA_peaks(1).peaks{k} ' subtracted'];
    header.name = filename; 
    save([filename '.mat'], 'data');
    save([filename '.lw6'], 'header');    
    clear data    
end
clear k e i data_i statement filename span

% choose EOIs based on centered datasets
for k = 1:length(GABA_peaks(1).peaks)
    % choose data, baseline correct
    data =  GABA_tracking(k).centered;
    
    % choose peak polarity
    if strcmp(GABA_peaks(1).peaks{k}(1), 'P')
        operation = 'max';
    elseif strcmp(GABA_peaks(1).peaks{k}(1), 'N')
        operation = 'min';
    end
    
    % loop through electrodes
    for e = 1:size(data, 1)
        % calculate maximal amplitude
        statement = ['peak_value(e) = ' operation '(data(e, :));'];
        eval(statement)
    end
    
    % identify three electrodes with the biggest response
    for e = 1:eoi_n
        statement = ['eoi_value(e) = ' operation '(peak_value);'];
        eval(statement)
        GABA_tracking(k).eoi.number(e) = find(peak_value == eoi_value(e));
        GABA_tracking(k).eoi.name{e} = labels(find(peak_value == eoi_value(e))); 
        peak_value(find(peak_value == eoi_value(e))) = 0;
    end
    
    % visualization parameters
    x_lim = ((1 + buffer) * GABA_peaks(1).widths(k))/2;
    x_eoi = -x_lim:xstep:x_lim;
    col = [0.8, 0.11, 0.23; 0.9, 0.24, 0.24; 1, 0.45, 0.45];
    
    % plot 
    fig = figure(figure_counter)
    hold on
    counter = 1;
    for e = 1:size(data, 1)
        if any(GABA_tracking(k).eoi.number == e)
            P(counter) = plot(x_eoi, data(e, :), 'Color', col(counter, :), 'LineWidth', 2.5);
            counter = counter + 1;
        else
            plot(x_eoi, data(e, :), 'Color', [0.65, 0.65, 0.65], 'LineWidth', 1.5)
        end
    end
    
    % add other parameters
    title(GABA_peaks(1).peaks{k})
    xlabel('time (s)')
    ylabel('amplitude (\muV)')
    set(gca, 'FontSize', 14)
    lgd = legend(P, {cell2mat(GABA_tracking(k).eoi.name{1}) cell2mat(GABA_tracking(k).eoi.name{2}) cell2mat(GABA_tracking(k).eoi.name{3})}, ...
        'Location', 'northeast');
    lgd.FontSize = 14;
    hold off
    
    % name and save figure
    figure_name = ['TEP_' target '_eoi'];
    savefig([figure_name])
    saveas(fig, [figure_name '.png'])
        
    % update figure counter
    figure_counter = figure_counter + 1 ;
end
clear k e operation statement fig x_lim x col figure_name lgd P 

% append new variables to the general MATLAB file
save(filename, 'GABA_tracking', '-append');
clear seed_electrode seed_peaks buffer eoi_n
                    
%% 5) peak windows
%% 6) amplitude extraction
%% functions
function peak_x = gmfp_plot(x, y, time_window, xstep, labeled, varargin)
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
[pks, locs] = findpeaks(y, 'MinPeakDistance', 10, 'MinPeakProminence', 0.015);
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

function [peak_x, peak_y] = gmfp_peaks(y, time_window, xstep, varargin)
% check whether to plot labels (default)
if ~isempty(varargin)
    a = find(strcmpi(varargin, 'max_peaks'));
    if ~isempty(a)
        max_peaks = varargin{a + 1};
    end
end

% find peaks 
[pks, locs] = findpeaks(y, 'MinPeakDistance', 10, 'MinPeakProminence', 0.015);
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
end

function [p1, p2, f] = gmfp_plot_diff(x, y, time_window, colours)
% calculate the time difference
for g = 1:size(y, 2)
    y(3, g) = y(2, g) - y(1, g);
end

% launch the figure 
hold on
plot(x, y(1, :))
plot(x, y(3, :))
yl = get(gca, 'ylim'); yl = [yl(1) - 0.1 , yl(2) + 1];
clf

% plot interpolated part
hold on
xlim(time_window)
rectangle('Position', [0, yl(1), 0.01, yl(2) - yl(1)], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')

% plot data
p1 = plot(x, y(1, :), 'Color', [0 0 0], 'LineWidth', 2.5);
p2 = plot(x, y(2, :), 'Color', colours(2, :) , 'LineWidth', 2.5);
f = fill([x fliplr(x)],[y(3, :) fliplr(zeros(1, length(x)))], colours(1, :), 'linestyle', 'none');
    
% add lines    
line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)
line(time_window, [0, 0], 'Color', [0, 0, 0], 'LineWidth', 1.5)

% add parameters
set(gca, 'fontsize', 14)
ylim(yl)
xlabel('time (s)')
ylabel('power (\muV^2)')
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

function [pos_x, data, sub_data] = track_peak(data, header, time_window, k, TEPs, buffer, seed)
% figure params 
figure_name = ['Peak ' TEPs.peaks{k}] ;
figure_center = TEPs.latencies(k);
span = ((1 + buffer) * TEPs.widths(k));

% set the peak timepoint manually
finish = 0;
while finish == 0
    % identify the TOI of current peak                    
    x1 = ceil((figure_center - span/2 - time_window(1)) / header.xstep);
    x2 = ceil((figure_center + span/2 - time_window(1)) / header.xstep);
    x = (figure_center - span/2) : header.xstep : (figure_center + span/2);

    % prepare data and header for visualization
    data_visual = data(1, :, :, :, :, [x1 : x2]);
    header_visual = header;
    header_visual.datasize(6) = length(data_visual);  
    header_visual.xstart = figure_center - span/2;

    % launch the figure                    
    axesHandles = [];
    fig = figure;   
    hold on          
    
    % check the length
    if length(squeeze(data_visual(1, seed, :, :, :, :))) > numel(x)
        delta = length(squeeze(data_visual(1, seed, :, :, :, :))) -  numel(x);
        data_visual = data_visual(1, :, :, :, :, 1:end-delta);
    elseif length(squeeze(data_visual(1, seed, :, :, :, :))) < numel(x)
        delta = numel(x) - length(squeeze(data_visual(1, seed, :, :, :, :)));
        data_visual(1, :, :, :, :, end + 1:end + delta) = zeros(1, delta);
    end
    
    % plot windowed data
    axesHandles = [axesHandles subplot(3, 3, [4:9])];  
    plot(x, squeeze(data_visual(1, seed, :, :, :, :)), 'LineWidth', 2)
    xlim([(figure_center - span/2), (figure_center + span/2)])
    title(figure_name, 'FontSize', 16)

    % plot the line at the center
    h = line([figure_center, figure_center], get(gca,'ylim'), 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--');   

    % plot the central topography 
    map_lims = get(axesHandles(1), 'ylim');
    axesHandles = [axesHandles subplot(3, 3, 2)];
    TEP_topoplot(header_visual, data_visual, figure_center, map_lims);

    % make the topography change with mouse movement 
    set(fig, 'WindowButtonMotionFcn', {@mouse_move, axesHandles, header_visual, data_visual});              

    % choose the peak position
    pos_x = get_position(axesHandles);             

    % update the figure
    set (fig, 'WindowButtonMotionFcn', '');
    subplot(3, 3, [4:9])
    set(h, 'XData', [pos_x, pos_x], 'YData', get(gca,'ylim'), 'LineStyle', '-');
    subplot(3, 3, 2) 
    cla(axesHandles(2))
    TEP_topoplot(header_visual, data_visual, pos_x, map_lims);

    pause(2)

    % ask for approval
    answer = questdlg('Do you want to proceed?', figure_name,...
        'Yes', 'No, I want to keep fiddling', 'Yes');

    % switch action
    switch answer
        case 'Yes'
            % close the figure
            close(fig)

            % exit the while loop
            finish = 1;

        case 'No, I want to keep fiddling'
            % assign previous center
            figure_center = pos_x;  

            % close the figure
            close(fig)                                
    end
end   

% index the peak timepoint  
n_peak = ceil((pos_x - time_window(1)) / header.xstep); 

% define final length
window_length = numel((pos_x - span/2) : header.xstep : (pos_x + span/2));

% loop through the data and subtract the value at the peak latency from all
% timepoints
sub_data = data;
topo_vector = squeeze(data(1, :, 1, 1, 1, n_peak));                                             
for i = 1:header.datasize(2)    % all channels
    for n = 1:length(data)      % all datapoints 
        sub_data(1, i, :, :, :, n) = data(1, i, :, :, :, n) - topo_vector(i);
    end
end

% crop the data
x_out1 = ceil((pos_x - span/2 - time_window(1)) / header.xstep);
x_out2 = ceil((pos_x + span/2 - time_window(1)) / header.xstep);
sub_data = squeeze(sub_data(1, :, :, :, :, [x_out1 : x_out2]));
data = squeeze(data(1, :, :, :, :, [x_out1 : x_out2]));

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

topoplot(vector2,chanlocs2,varargin{:});
set(gcf,'color',[1 1 1]);
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
