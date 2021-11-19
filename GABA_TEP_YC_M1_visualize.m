%% GABA-AD: TMS-EVOKED POTENTIALS - PRIMARY MOTOR CORTEX
% Written by Dominika for GABA-AD project (2020-21)
% 
% Colection of scripts to visualy explore extracted TMS-EEG variables: 
%   --> output datasets are saved in a folder 'GABA_YC_statistics'
%   --> figures are saved in a folder 'GABA_YC_figures'
% 
% 1) load the data


%% parameters
clear all; clc

% ----- adjustable parameters -----
% dataset
prefix = 'GABA';
group = 'YC';
participant = 1:20;
medication = {'placebo' 'alprazolam'}; 
time = {'pre' 'post'};
stimulus = {'CS' 'TS' 'ppTMS'};

% statistics
z = 1.96;
alpha = 0.05;
% --------------------------------

% navigate to the Git folder --> source of saved default files
folder_git = uigetdir(pwd, 'Choose the Git folder');

% visualization 
figure_counter = 1;

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
folder_results = uigetdir(pwd, 'Choose the Results folder');
folder_input = [folder_results '\GABA_' group '_variables'];
folder_output = [folder_results '\GABA_' group '_statistics'];
folder_figures = [folder_results '\GABA_' group '_figures'];
if ~exist(folder_output)
    mkdir(folder_output)
end
if ~exist(folder_figures)
    mkdir(folder_figures)  
end
clear folder_results

%% 1) load data
% TEPs from M1
load([folder_input '\' prefix '_' group '_M1_TEPs.mat'], 'GABA_TEP_peaks', 'GABA_TEP_default');

%% 2) plot histograms per TEP component
% peak amplitude
for p = 1:length(GABA_TEP_default.peak)
    % launch the figure
    fig = figure(figure_counter);
    
    % baseline
    for m = 1:length(medication)
        for s = 1:length(stimulus)
            % plot the histogram
            subplot(length(medication)*length(time), length(stimulus), (m-1)*3 + s)
            histogram(GABA_TEP_peaks.amplitude_peak(m, 1, s, :, p), 10, 'FaceColor', colours((m-1)*2 + 1, :), 'EdgeColor', [1 1 1])
            
            % set annotation
            if s == 1
                yl = get(gca, 'ylim');
                xl = get(gca, 'xlim');
                text(xl(1)- (xl(2) - xl(1))*0.5, yl(2)/2 + yl(2)*0.1, sprintf('baseline\n%s', medication{m}), 'FontSize', 10)
                clear yl xl
            end
                       
            % set title
            if m == 1
                title(stimulus{s})
            end
        end 
    end
    
    % post medication
    for m = 1:length(medication)
        for s = 1:length(stimulus)
            % plot the histogram
            subplot(length(medication)*length(time), length(stimulus), (m+1)*3 + s)
            histogram(GABA_TEP_peaks.amplitude_peak(m, 2, s, :, p), 10, 'FaceColor', colours((m-1)*2 + 2, :), 'EdgeColor', [1 1 1])
            
            % set annotation
            if s == 1
                yl = get(gca, 'ylim');
                xl = get(gca, 'xlim');
                text(xl(1)- (xl(2) - xl(1))*0.5, yl(2)/2 + yl(2)*0.1, sprintf('post\n%s', medication{m}), 'FontSize', 10)
                clear yl xl
            end            
        end 
    end
    
    % parameters
    sgtitle([GABA_TEP_default.peak{p} ' - raw data'])
        
    % update counter
    figure_counter = figure_counter + 1;
end
clear p fig m s

% define skewness
skewness = [-1, 1, 1, 1, -1, 1];

%% 3) square root
% add a constant to make all entries positive
c = floor(min(GABA_TEP_peaks.amplitude_peak, [], 'all'));
GABA_TEP_peaks.amplitude_peak_sqrt = GABA_TEP_peaks.amplitude_peak - c;

% reflect negatively skewed data
for p = find(skewness == -1)
    for m = 1:length(medication)
        for t = 1:length(time)
            for s = 1:length(stimulus)
                GABA_TEP_peaks.amplitude_peak_sqrt(m, t, s, :, p) = max(GABA_TEP_peaks.amplitude_peak_sqrt(m, t, s, :, p)) - GABA_TEP_peaks.amplitude_peak_sqrt(m, t, s, :, p) + 1;
            end
        end
    end
end
clear p m t s

% square-root the data
GABA_TEP_peaks.amplitude_peak_sqrt = sqrt(GABA_TEP_peaks.amplitude_peak_sqrt);

% replot
for p = 1:length(GABA_TEP_default.peak)
    % launch the figure
    fig = figure(figure_counter);
    
    % baseline
    for m = 1:length(medication)
        for s = 1:length(stimulus)
            % plot the histogram
            subplot(length(medication)*length(time), length(stimulus), (m-1)*3 + s)
            histogram(GABA_TEP_peaks.amplitude_peak_sqrt(m, 1, s, :, p), 10, 'FaceColor', colours((m-1)*2 + 1, :), 'EdgeColor', [1 1 1])
            
            % set annotation
            if s == 1
                yl = get(gca, 'ylim');
                xl = get(gca, 'xlim');
                text(xl(1)- (xl(2) - xl(1))*0.5, yl(2)/2 + yl(2)*0.1, sprintf('baseline\n%s', medication{m}), 'FontSize', 10)
                clear yl xl
            end
                       
            % set title
            if m == 1
                title(stimulus{s})
            end
        end 
    end
    
    % post medication
    for m = 1:length(medication)
        for s = 1:length(stimulus)
            % plot the histogram
            subplot(length(medication)*length(time), length(stimulus), (m+1)*3 + s)
            histogram(GABA_TEP_peaks.amplitude_peak_sqrt(m, 2, s, :, p), 10, 'FaceColor', colours((m-1)*2 + 2, :), 'EdgeColor', [1 1 1])
            
            % set annotation
            if s == 1
                yl = get(gca, 'ylim');
                xl = get(gca, 'xlim');
                text(xl(1)- (xl(2) - xl(1))*0.5, yl(2)/2 + yl(2)*0.1, sprintf('post\n%s', medication{m}), 'FontSize', 10)
                clear yl xl
            end            
        end 
    end
    
    % parameters
    sgtitle([GABA_TEP_default.peak{p} ' - square root of the data'])
        
    % update counter
    figure_counter = figure_counter + 1;
end
clear p fig m s

% save data in a R-compatible table 
GABA_TEP_peaks_sqrt = table;
row_counter = 1;
for p = 1:length(participant) 
    for m = 1:length(medication)  
        for t = 1:length(time)
            for s = 1:length(stimulus)
                for k = 1:length(GABA_TEP_default.peak) 
                    %fill in the table
                    GABA_TEP_peaks_sqrt.subject(row_counter) = participant(p);
                    GABA_TEP_peaks_sqrt.medication(row_counter) = medication(m);
                    GABA_TEP_peaks_sqrt.time(row_counter) = time(t);
                    GABA_TEP_peaks_sqrt.stimulus(row_counter) = stimulus(s);
                    GABA_TEP_peaks_sqrt.peak(row_counter) = GABA_TEP_default.peak(k);
                    GABA_TEP_peaks_sqrt.amplitude_peak(row_counter) = GABA_TEP_peaks.amplitude_peak_sqrt(m, t, s, p, k);
                    GABA_TEP_peaks_sqrt.latency(row_counter) = GABA_TEP_peaks.latency(m, t, s, p, k);
                    
                    % update the counter
                    row_counter = row_counter + 1;
                end
            end
        end
    end
end
clear m t s p k row_counter
writetable(GABA_TEP_peaks_sqrt, [folder_output '\GABA_TEP_peaks_sqrt.csv'])

%% 4) log-transform
% add a constant to make all entries positive
c = floor(min(GABA_TEP_peaks.amplitude_peak, [], 'all'));
GABA_TEP_peaks.amplitude_peak_log = GABA_TEP_peaks.amplitude_peak - c;

% reflect negatively skewed data
for p = find(skewness == -1)
    for m = 1:length(medication)
        for t = 1:length(time)
            for s = 1:length(stimulus)
                GABA_TEP_peaks.amplitude_peak_log(m, t, s, :, p) = max(GABA_TEP_peaks.amplitude_peak_log(m, t, s, :, p)) - GABA_TEP_peaks.amplitude_peak_log(m, t, s, :, p) + 1;
            end
        end
    end
end
clear p m t s

% slog transform the data
GABA_TEP_peaks.amplitude_peak_log = log10(GABA_TEP_peaks.amplitude_peak_log);

% replot
for p = 1:length(GABA_TEP_default.peak)
    % launch the figure
    fig = figure(figure_counter);
    
    % baseline
    for m = 1:length(medication)
        for s = 1:length(stimulus)
            % plot the histogram
            subplot(length(medication)*length(time), length(stimulus), (m-1)*3 + s)
            histogram(GABA_TEP_peaks.amplitude_peak_log(m, 1, s, :, p), 10, 'FaceColor', colours((m-1)*2 + 1, :), 'EdgeColor', [1 1 1])
            
            % set annotation
            if s == 1
                yl = get(gca, 'ylim');
                xl = get(gca, 'xlim');
                text(xl(1)- (xl(2) - xl(1))*0.5, yl(2)/2 + yl(2)*0.1, sprintf('baseline\n%s', medication{m}), 'FontSize', 10)
                clear yl xl
            end
                       
            % set title
            if m == 1
                title(stimulus{s})
            end
        end 
    end
    
    % post medication
    for m = 1:length(medication)
        for s = 1:length(stimulus)
            % plot the histogram
            subplot(length(medication)*length(time), length(stimulus), (m+1)*3 + s)
            histogram(GABA_TEP_peaks.amplitude_peak_log(m, 2, s, :, p), 10, 'FaceColor', colours((m-1)*2 + 2, :), 'EdgeColor', [1 1 1])
            
            % set annotation
            if s == 1
                yl = get(gca, 'ylim');
                xl = get(gca, 'xlim');
                text(xl(1)- (xl(2) - xl(1))*0.5, yl(2)/2 + yl(2)*0.1, sprintf('post\n%s', medication{m}), 'FontSize', 10)
                clear yl xl
            end            
        end 
    end
    
    % parameters
    sgtitle([GABA_TEP_default.peak{p} ' - log10 of the data'])
        
    % update counter
    figure_counter = figure_counter + 1;
end
clear p fig m s

% save data in a R-compatible table 
GABA_TEP_peaks_log = table;
row_counter = 1;
for p = 1:length(participant) 
    for m = 1:length(medication)  
        for t = 1:length(time)
            for s = 1:length(stimulus)
                for k = 1:length(GABA_TEP_default.peak) 
                    %fill in the table
                    GABA_TEP_peaks_log.subject(row_counter) = participant(p);
                    GABA_TEP_peaks_log.medication(row_counter) = medication(m);
                    GABA_TEP_peaks_log.time(row_counter) = time(t);
                    GABA_TEP_peaks_log.stimulus(row_counter) = stimulus(s);
                    GABA_TEP_peaks_log.peak(row_counter) = GABA_TEP_default.peak(k);
                    GABA_TEP_peaks_log.amplitude_peak(row_counter) = GABA_TEP_peaks.amplitude_peak_log(m, t, s, p, k);
                    GABA_TEP_peaks_log.latency(row_counter) = GABA_TEP_peaks.latency(m, t, s, p, k);
                    
                    % update the counter
                    row_counter = row_counter + 1;
                end
            end
        end
    end
end
clear m t s p k row_counter
writetable(GABA_TEP_peaks_log, [folder_output '\GABA_TEP_peaks_log.csv'])