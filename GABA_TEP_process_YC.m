%% GABA-AD: TMS-EVOKED POTENTIALS
% Written by Dominika for GABA-AD project (2020-21)
% 
% Colection of scripts to analyse preprocessed TMS-EEG data: 
% 1) calculates global mean field power from grand average
%       - averages baseline and post-medication data from both sessions
%       - calculates GMFP and plots it, saves the figure
%       - automatically identifies local peaks within [0.01 0.25]s and
%       saves peak latencies in vector 'TEP_peaks'
% 2) determines electrodes of interest
%       --> average across all conditions for each subjects
%       - uses GMFP peak latencies as default search points and signal from Cz 
%       as starting peak point to track peak progress
%       - averages across subjects
%       - identifies 3 electrodes with the highest average amplitude for
%       each peak 
%       - pools identified electrodes into an EOI for each peak 
%       - saves new individual datasets with EOIs for all analysed peaks
% 3) determines peak windows 
%       --> average across all conditions and all subjects, EOI signals
%       - calculates half prominence of each analysed peak
% 4) extracts mean amplitude  


%% parameters
clear all; clc

% dataset
group = 'YC';
participant = 1:20;
medication = {'placebo' 'alprazolam'}; 
time = {'pre' 'post'};

% target --> stimulus, labels, prefix, header
target = questdlg('Choose target:', 'Stimulated cortex', 'M1', 'AG', 'M1');
load('labels.mat')
switch target
    case 'M1'
        stimulus = {'CS' 'TS' 'ppTMS'};
        labels.CS = labels_CS; labels.TS = labels_TS;
        prefix = 'eois crop avg final_dataset';
        load([prefix ' ' stimulus{1} ' ' group num2str(participant(1)) ' ' medication{1} ' ' time{1} ' ' target '.lw6'], '-mat')
    case 'AG'
        stimulus = {''};
        labels = labels_AG;     
        prefix = 'avg bl icfilt ica26 crop but fft-notchfilt prefilt prea28 visual reref ds art-sup bl ep dc chan-select';
        load([prefix ' ' group num2str(participant(1)) ' ' medication{1} ' ' target ' ' time{1} '.lw6'], '-mat')
end
clear labels_CS labels_TS labels_AG

% visualization inserted params 
time_window = [-0.05, 0.3];
z = 1.96;
alpha = 0.2;

% visualization calculated params
figure_counter = 1;
xstep = header.xstep; 
x = [time_window(1):xstep:time_window(2)];
x_start = (time_window(1) - header.xstart)/xstep;
x_end = (time_window(2) - header.xstart)/xstep;
load('colours.mat'); colours = colours([2 1 4 6],:);
clear header

%% 1) GMFP
% decide output parameters 
labeled = 'on';
max_peaks = 5;

% all conditions together
for s = 1:length(stimulus)
    % load grand average data
    a = 1;
    for m = 1:length(medication)
        for t = 1:length(time)
            switch target
                case 'AG'
                    load(['avg merged ' medication{m} ' ' time{t} ' ' target '.mat'])
                case 'M1'
                    load(['avg merged ' medication{m} ' ' time{t} ' ' stimulus{s} ' ' target '.mat'])
            end
            gmfp_data(a, :, :) = squeeze(data(:, 1:30, :, :, :, x_start:x_end));
            a = a + 1;
        end 
    end
    clear a m t data 

    % average conditions 
    for a = 1:size(gmfp_data, 2)
        for b = 1:size(gmfp_data, 3)
            gmfp(a, b) = mean(gmfp_data(:, a, b));
        end
    end
    clear a b gmfp_data
    
    % dataset name + figure title
    switch target
        case 'AG'
            name = ['GMFP_' target '_merged'];
            fig_title = ['GMFP : ' target];
        case 'M1'
            name = ['GMFP_' target '_' stimulus{s} '_merged'];
            fig_title = ['GMFP : ' target ', ' stimulus{s}];
    end

    % calculate GMFP
    gmfp = std(gmfp, 1);  
    GMFP(s, :) = gmfp;
    
    % plot GMFP and extract peak latencies
    fig = figure(figure_counter);
    if ~isempty(max_peaks)
        TEP_peaks(s, :) = plot_gmfp(x, gmfp, time_window, xstep, labeled, 'max_peaks', max_peaks);
    else
        TEP_peaks(s, :) = plot_gmfp(x, gmfp, time_window, xstep, labeled);
    end
    title(fig_title, 'fontsize', 16, 'fontweight', 'bold')

    % save figure, update    
    savefig(name)
    saveas(fig, [name '.png'])
    figure_counter = figure_counter + 1;
end
save(['GMFP_' target '_merged.mat'], 'GMFP')
clear s GMFP gmfp gmfp_data fig name 

% medication and time facors separately
% s = 1; m = 1; t = 1;
for s = 1:length(stimulus)
    for m = 1:length(medication)
        for t = 1:length(time)
            % load grand average data                
            switch target
                case 'AG'
                    load(['avg merged ' medication{m} ' ' time{t} ' ' target '.mat'])
                case 'M1'
                    load(['avg merged ' medication{m} ' ' time{t} ' ' stimulus{s} ' ' target '.mat'])
            end
            gmfp_data = squeeze(data(:, 1:30, :, :, :, x_start:x_end));
            clear data 

            % calculate GMFP
            gmfp(s, m, t, :) = std(gmfp_data, 1);
            clear gmfp_data 
        end
        
        % dataset name + figure title
        switch target
            case 'AG'
                name = ['GMFP_' target '_' medication{m}];
                fig_title = ['GMPF - ' target ' : ' medication{m} ];
            case 'M1'
                name = ['GMFP_' target '_' medication{m} '_' stimulus{s}];
                fig_title = ['GMPF - ' target ', ' stimulus{s} ' : ' medication{m} ];
        end
        
        % plot GMFP per conditions (baseline X post medication)
        fig = figure(figure_counter);
        [p1, p2, f] = plot_gmfp_diff(x, squeeze(gmfp(s, m, :, :)), time_window, colours([(m - 1)*2 + 1, (m - 1)*2 + 2], :));
        title(fig_title, 'fontsize', 16, 'fontweight', 'bold')
        legend([p1, p2, f], {'baseline' 'post medication' 'change'}, 'FontSize', 12, ...
            'Location', 'northeast', 'NumColumns', 2, 'EdgeColor', 'none');

        % save figure, update    
        savefig(name)
        saveas(fig, [name '.png'])
        figure_counter = figure_counter + 1;
    end
end
save(['GMFP_' target '_conditions.mat'], 'gmfp')
clear s m t gmfp fig name p1 p2 f  

% functions
function peak_x = plot_gmfp(x, y, time_window, xstep, labeled, varargin)
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
clf

% plot interpolated part
hold on
xlim(time_window)
rectangle('Position', [0, yl(1), 0.01, yl(2) - yl(1)], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')

% plot data, mark TMS stimulus
plot(x, y, 'Color', [0 0 0], 'LineWidth', 2.5)
line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)

% find peaks 
[pks, locs] = findpeaks(y);
for a = 1:length(locs)
    if time_window(1) + locs(a)*xstep <= 0.01
        idx(a) = false;
    elseif time_window(1) + locs(a)*xstep > 0.25
        idx(a) = false;        
    else
        idx(a) = true;
    end
end
pks = pks(idx); locs = locs(idx);
pks = pks(1:max_peaks); locs = locs(1:max_peaks);

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
ylabel('power (\muV^2)')
end
function [p1, p2, f] = plot_gmfp_diff(x, y, time_window, colours)
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
%% 2) EOIs
%% 3) peak windows
%% 4) amplitude


