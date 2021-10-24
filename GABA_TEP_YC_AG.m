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

% dataset parameters
load([prefix ' ' group ' 01 ' target ' ' medication{1} ' ' time{1} '.lw6'], '-mat')
labels = {header.chanlocs.labels};

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
folderpath = uigetdir;
if ~exist([folderpath '\GABA_' group '_export'])
    mkdir(folderpath, ['GABA_' group '_export'])
elseif ~exist([folderpath '\GABA_' group '_figures'])
    mkdir(folderpath, ['GABA_' group '_figures'])  
end

%% 1) extract individual data
% load data --> uncorrected ppTMS TEPs
for m = 1:length(medication)
    for t = 1:length(time)
        for p = 1:length(participant)
            % define participant
            if p < 10
                subj = ['0' num2str(participant(p))];
            else
                subj = num2str(participant(p));
            end

            % load individual dataset
            load([prefix ' ' group ' ' subj ' ' target ' ' medication{m} ' ' time{t} '.mat'])

            % append the data in the data matrix
            GABA_data(m, t, p, :, :) = squeeze(data(:, :, :, :, :, x_start:x_end));                                  
        end
    end
end
clear m t p data subj
disp(['Data import finished. Datasize: ' num2str(size(GABA_data))])

% save dataset to the global MATLAB file
if ~exist([folderpath '\' filename '.mat']) 
    save([folderpath '\' filename '.mat'], 'GABA_data');
else
    save([folderpath '\' filename '.mat'], 'GABA_data', '-append');
end

%% 2) preliminary visualization 
% ----- decide output parameters -----
electrode = {'Cz'};
% ------------------------------------
% load data if not loaded
if exist('GABA_data') ~= 1
    load([folderpath '\' filename '.mat'], 'GABA_data')
end

% prepare average statistics
for m = 1:length(medication)
    for t = 1:length(time)
        for e = 1:size(GABA_data, 4)
            for i = 1:size(GABA_data, 5)
                data_mean(m, t, e, i) = mean(squeeze(GABA_data(m, t, :, e, i)));
                data_CI(m, t, e, i) = (std(squeeze(GABA_data(m, t, :, e, i)))/sqrt(length(participant))) * z;
            end
        end
    end
end
clear m t e i

% plot the baseline data both sessions
for e = 1:length(electrode) 
    % identify the electrode
    e_n = find(contains(labels, electrode{e}));
    
    % prepare data
    data_visual = squeeze(data_mean(:, 1, e_n, :));
    CI_visual = squeeze(data_CI(:, 1, e_n, :));

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
    title([target ', ' electrode{e} ' electrode: baseline TEP'])
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
    savefig([folderpath '\GABA_' group '_figures\' figure_name '.fig'])
    saveas(fig, [folderpath '\GABA_' group '_figures\' figure_name '.png'])

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
        data_visual = squeeze(data_mean(m, :, e_n, :));
        CI_visual = squeeze(data_CI(m, :, e_n, :));

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
        title([target ', ' electrode{e} ' electrode: ' medication{m}])
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
        savefig([folderpath '\GABA_' group '_figures\' figure_name '.fig'])
        saveas(fig, [folderpath '\GABA_' group '_figures\' figure_name '.png'])

        % update figure counter
        figure_counter = figure_counter + 1 ;
    end
end
clear m e t e_n fig lgd data_visual CI_visual figure_name P F yl

% save dataset to the global MATLAB file
save([folderpath '\' filename '.mat'], 'data_mean', 'data_CI', '-append');
clear electrode 

%% 3) GFP - peak identification
% ----- decide output parameters -----
labeled = 'off';
max_peaks = 6;
% ------------------------------------
if exist('GABA_data') ~= 1
    load([folderpath '\' filename '.mat'], 'GABA_data')
end

% compute individual GFP
GABA_GFP_subject = struct;
for m = 1:length(medication)
    for t = 1:length(time)
        for p = 1:length(participant)
            % calculate GFP (exclude target channel)
            GABA_GFP_subject.data(m, t, p, :) = std(squeeze(GABA_data(m, t, p, 1:30, :)), 1);  
        end
    end
end
clear m t p 

% compute grand average GFP, plot
row_count = 1;
for m = 1:length(medication)
    for t = 1:length(time)
        % calculate GFP 
        GABA_GFP_mean(m, t, :) = std(squeeze(data_mean(m, t, 1:30, :)), 1);  

        % set dataset name + figure title
        fig_name = ['GABA_GFP_' medication{m} '_' time{t}];
        fig_title = ['GFP: ' medication{m} ', ' time{t} ' medication'];

        % plot GFP and extract peak latencies
        fig = figure(figure_counter);
        hold on

        if ~isempty(max_peaks)
        % plot GFP
        h_axis(1) = subplot(3, max_peaks, [1 : 2*max_peaks]);
        GABA_TEP_default(row_count).latencies = gfp_plot(x, squeeze(GABA_GFP_mean(m, t, :)), time_window, xstep, labeled, 'max_peaks', max_peaks);
        title(fig_title, 'fontsize', 16, 'fontweight', 'bold')

            % add topoplots
            for t = 1:length(GABA_TEP_default(row_count).latencies)
                % choose data for topoplot 
                for e = 1:size(GABA_data, 4)
                    data_topoplot(1, e, 1, 1, 1, :) = squeeze(gfp_data(p, c, e, :));
                end

                % plot the topoplot
                h_axis(1 + t) = subplot(3, max_peaks, 2*max_peaks + t);            
                topo_plot(header, data_topoplot, GABA_TEP_default(row_count).latencies(t), time_window(1), [-2, 2])

                % shift down
                pos = get(h_axis(1 + t), 'Position');
                pos(2) = pos(2) - 0.05;
                set(h_axis(1 + t), 'Position', pos);

                % add timing
                text(-0.3, -0.8, sprintf('%1.0f ms', AGSICI_TEP(row_count).latencies(t)*1000), 'Color', [1 0 0], 'FontSize', 14)
            end

        else
            GABA_TEP_default(row_count).latencies(:) = gfp_plot(x, squeeze(AGSICI_GFP(p, c, :)), time_window, xstep, labeled);
        end  
        hold off

        % save figure, update    
        savefig([pwd '\' filename '_figs\' fig_name '.fig'])
        saveas(fig, [pwd '\' filename '_figs\' fig_name '.png'])
        figure_counter = figure_counter + 1;
        
        % update row counter
        row_count = row_count + 1;
    end
end
clear p c t h_axis pos fig_name fig_title row_count fig data_topoplot

% append new variables to the general MATLAB file
save([filename '.mat'], 'AGSICI_GFP', 'AGSICI_TEP', 'AGSICI_GFP_subject', '-append')
clear labeled max_peaks gfp_CI

% append new variables to the general MATLAB file
save([folderpath '\' filename '.mat'], 'GABA_GFP_subject', 'GABA_GFP_mean', '-append');
clear labeled max_peaks 

%% 4) GFP - plot difference
% ----- decide output parameters -----
TOI_peaks = [0.03 0.045 0.060 0.100 0.180];
peaks = {'P30' 'N45' 'P60' 'N100' 'P180'};
% ------------------------------------
% calculate GFP for each subject/condition
for m = 1:length(medication)
    for t = 1:length(time)
        for s = 1:length(stimulus) 
            for p = 1:length(participant)
                % calculate GMFP
                GABA_GFP_subject(m, t, s, p, :) = std(squeeze(GABA_data(m, t, s, p, 1:30, :)), 1);
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
        data_visual = squeeze(mean(squeeze(GABA_GFP_subject(m, :, s, :, :)), 2));

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
save([folderpath '\' filename '.mat'], 'GABA_GFP_subject', '-append');
clear TOI TOI_peaks peaks

%% 5) DISS
% load data if not loaded
if exist('GABA_GFP_subject') ~= 1
    load([folderpath '\' filename '.mat'], 'GABA_GFP_subject')
end

% normalize data by GFP
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
% write text files for Ragu --> uncorrected ppTMS TEPs
for m = 1:length(medication)
    for t = 1:length(time)
        for p = 1:length(participant)
            % choose data to write, remove 'eoi' channels
            data = squeeze(GABA_data(m, t, p, 1:30, :))';

            % define subject
            if participant(p) < 10
                subj = ['S0' num2str(participant(p))];
            else
                subj = ['S' num2str(participant(p))];
            end

            % save as .csv               
            name = ['GABA_' group '_' subj '_' target '_' medication{m} '_' time{t} '.csv']; 
            writematrix(data, [folderpath '\GABA_' group '_export\' name])
        end
    end
end
clear m t p data subj name     

% create the montage file
name = [folderpath '\GABA_' group '_export\GABA_montage.xyz'];
fileID = fopen(name, 'a');
fprintf(fileID, '30\r\n');
for a = 1:30
    fprintf(fileID, '%.4f %.4f %.4f %s\r\n', ...
        header.chanlocs(a).X, header.chanlocs(a).Y, header.chanlocs(a).Z, header.chanlocs(a).labels);
end
fclose(fileID)
clear name fileID a folder_name
%% 6) amplitude extraction
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
ylabel('GFP (\muV)')
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