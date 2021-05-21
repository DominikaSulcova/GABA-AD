%% GABA-AD: TMS-EVOKED POTENTIALS
% Written by Dominika for GABA-AD project (2020-21)
% 
% % Colection of scripts to analyse preprocessed TMS-EEG data: 
% 1) load the data
%       - loads individual averaged data and trims it in a predefined time
%       window
%       - saves the data from all subjects and conditions in a 6D matrix:
%       2 (medication) X 2 (time) X number of stimuli X number of subjects X
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
participant = 1:20;
medication = {'placebo' 'alprazolam'}; 
time = {'pre' 'post'};

% visualization 
time_window = [-0.05, 0.3];
z = 1.96;
alpha = 0.2;

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
        prefix = 'eoi avg bl icfilt ica26 crop but fft-notchfilt prefilt prea28 visual reref ds art-sup bl ep dc chan-select';
        load([prefix ' ' group num2str(participant(1)) ' ' medication{1} ' ' target ' ' time{1} '.lw6'], '-mat')
end
clear labels_CS labels_TS labels_AG
% --------------------------------

% visualization calculated params
figure_counter = 1;
xstep = header.xstep; 
x = [time_window(1):xstep:time_window(2)];
x_start = (time_window(1) - header.xstart)/xstep;
x_end = (time_window(2) - header.xstart)/xstep;
load('colours.mat'); colours = colours([2 1 4 6],:);

%% 1) extract individual data
% load data based on the prefix + conditions
switch target
    case 'M1'
        for m = 1:length(medication)
            for t = 1:length(time)
                for s = 1:length(stimulus)
                    for p = 1:length(participant)
                        % load individual dataset
                        load([prefix ' ' stimulus{s} ' ' group num2str(participant(p)) ' ' medication{m} ' ' time{t} ' ' target '.mat'])
                        
                        % append the data in the data matrix
                        GABA_data(m, t, s, p, :, :) = squeeze(data(:, :, :, :, :, x_start:x_end));                        
                    end
                end
            end
        end
    case 'AG'
        for m = 1:length(medication)
            for t = 1:length(time)
                for s = 1:length(stimulus)
                    for p = 1:length(participant)
                        % load individual dataset
                        load([prefix ' ' group num2str(participant(p)) ' ' medication{m} ' ' target ' ' time{t} '.mat'])
                        
                        % append the data in the data matrix
                        GABA_data(m, t, s, p, :, :) = squeeze(data(:, :, :, :, :, x_start:x_end));                        
                    end
                end
            end
        end
end
clear m t s p data
disp(['Data import finished. Datasize: ' num2str(size(GABA_data))])

% save dataset to the global MATLAB file
filename = ['GABA_' target '.mat'];
if exist(filename) == 0
    save(filename, 'GABA_data');
else
    save(filename, 'GABA_data', '-append');
end

%% 2) preliminary visualization 
% ----- decide output parameters -----
electrode = {'target'};
% ------------------------------------

% prepare baseline statistics
for m = 1:length(medication)
    for e = 1:size(GABA_data, 5)
        for i = 1:size(GABA_data, 6)
            baseline_mean(m, e, i) = mean(squeeze(GABA_data(m, 1, 1, :, e, i)));
            baseline_CI(m, e, i) = (std(squeeze(GABA_data(m, 1, 1, :, e, i)))/sqrt(length(participant))) * z;
        end
    end
end
clear m e i

% plot the data
for e = 1:length(electrode)  
    % prepare data
    data_visual = squeeze(baseline_mean(:, find(contains(labels, electrode{e})), :));
    CI_visual = squeeze(baseline_CI(:, find(contains(labels, electrode{e})), :));
    
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
    lgd = legend(P, medication, 'Location', 'southwest');
    lgd.FontSize = 14;
    hold off
    
    % name and save figure
    figure_name = ['TEP_' target '_baseline_' electrode{e}];
    savefig([figure_name])
    saveas(fig, [figure_name '.png'])
        
    % update figure counter
    figure_counter = figure_counter + 1 ;
end
clear m e fig lgd data_visual CI_visual figure_name P F yl

% save dataset to the global MATLAB file
save(filename, 'baseline_mean', 'baseline_CI', '-append');
clear electrode baseline_mean baseline_CI 

%% 3) GMFP
% ----- decide output parameters -----
labeled = 'off';
max_peaks = 5;
peak_names = {'P20' 'N35' 'P50' 'P75' 'N100' 'P180'};     % choose names of peaks based on the visual inspection of average signal
% ------------------------------------

% all conditions together + topoplots
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
            gmfp(a, b) = mean(squeeze(gmfp_data(:, a, b)));
        end
    end
    clear a b gmfp_data
    
    % dataset name + figure title
    switch target
        case 'AG'
            name = ['GMFP_' target '_merged'];
            fig_title = 'GMFP : angular gyrus';
        case 'M1'
            name = ['GMFP_' target '_' stimulus{s} '_merged'];
            fig_title = ['GMFP : primary motor cortex, ' stimulus{s}];
    end

    % calculate GMFP
    GABA_GMFP(s, :) = std(gmfp, 1);  
    
    % plot GMFP and extract peak latencies
    fig = figure(figure_counter);
    hold on
    
    if ~isempty(max_peaks)
        % plot gmfp
        h_axis(1) = subplot(3, max_peaks, [1 : 2*max_peaks]);
        GABA_TEP(s).peaks = peak_names;
        GABA_TEP(s).latencies = gmfp_plot(x, GABA_GMFP(s, :), time_window, xstep, labeled, 'max_peaks', max_peaks);
        title(fig_title, 'fontsize', 16, 'fontweight', 'bold')
        
        % add topoplots
        for t = 1:length(GABA_TEP(s).peaks)
            % plot the topoplot
            h_axis(1 + t) = subplot(3, max_peaks, 2*max_peaks + t);
            gmfp_topoplot(header, gmfp, GABA_TEP(s).latencies(t), time_window(1), [-2, 2])
            
            % shift down
            pos = get(h_axis(1 + t), 'Position');
            pos(2) = pos(2) - 0.05;
            set(h_axis(1 + t), 'Position', pos);
            
            % add timing
            text(-0.3, -0.8, sprintf('%1.0f ms', GABA_TEP(s).latencies(t)*1000), 'Color', [1 0 0], 'FontSize', 14)
        end
        
    else
        GABA_TEP(s).latencies = gmfp_plot(x, GABA_GMFP(s, :), time_window, xstep, labeled);
    end  
    hold off

    % save figure, update    
    savefig(name)
    saveas(fig, [name '.png'])
    figure_counter = figure_counter + 1;
end
clear s GMFP gmfp gmfp_data fig name pos h_axis

%%% ONLY FOR AG %%%
GABA_TEP(s).latencies = [GABA_TEP(s).latencies(1:2), 0.05, GABA_TEP(s).latencies(3:5)];
%%%%%%%%%%%%%%%%%%%

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
        [p1, p2, f] = gmfp_plot_diff(x, squeeze(gmfp(s, m, :, :)), time_window, colours([(m - 1)*2 + 1, (m - 1)*2 + 2], :));
        title(fig_title, 'fontsize', 16, 'fontweight', 'bold')
        legend([p1, p2, f], {'baseline' 'post medication' 'change'}, 'FontSize', 12, ...
            'Location', 'northeast', 'NumColumns', 2, 'EdgeColor', 'none');

        % save figure, update    
        savefig(name)
        saveas(fig, [name '.png'])
        figure_counter = figure_counter + 1;
    end
end
GABA_GMFP_conditions = gmfp; 
clear s m t gmfp fig name p1 p2 f fig_title

% append new variables to the general MATLAB file
save(filename, 'GABA_GMFP', 'GABA_TEP', 'GABA_GMFP_conditions', '-append');
clear labeled max_peaks peak_names

%% 4) peak widths
% calculate narrow window parameters
x_start_narrow = (0.01 - time_window(1))/xstep;
x_end_narrow = (0.25 - time_window(1))/xstep;
x_narrow = [0.01:xstep:0.25];

% loop through stimuli
for s = 1:length(stimulus)
    % choose data and x
    data = GABA_GMFP(s, x_start_narrow:x_end_narrow);
        
    % identify peak widths
    [P, L, GABA_TEP(s).widths, R] = findpeaks(data, 'Annotate','extents', 'WidthReference', 'halfheight');
    GABA_TEP(s).widths = ceil(GABA_TEP(s).widths)* header.xstep;
    if length(GABA_TEP(s).widths) ~= length(GABA_TEP(s).peaks)
        disp('ATTENTION: number of identified peaks does not match with peaks extracted in the previous step!')
    end
    
    % plot figure
    fig = figure(figure_counter);
    hold on
    findpeaks(data, x_narrow, 'Annotate', 'extents', 'WidthReference', 'halfheight');
    set(gca, 'fontsize', 14)
    xlabel('time(s)')
    ylabel('power (\muV^2)')
    grid off
    
    % add width denotation
    for k = 1:length(GABA_TEP(s).widths)
        if k == length(GABA_TEP(s).widths)
            text(GABA_TEP(s).latencies(k) - 0.005, -0.25, sprintf('%1.0f ms', GABA_TEP(s).widths(k)*1000), 'Color', [0.93 0.69 0.13], 'FontSize', 14)
        else
            text(GABA_TEP(s).latencies(k) - 0.005, -0.25, sprintf('%1.0f', GABA_TEP(s).widths(k)*1000), 'Color', [0.93 0.69 0.13], 'FontSize', 14)
        end
    end
    xlim([-0.005, 0.26])
    
    % save figure, update
    savefig(['GABA_' target '_widths'])
    saveas(fig, ['GABA_' target '_widths.png'])
    figure_counter = figure_counter + 1;
end

%%% ONLY FOR AG %%%
GABA_TEP(s).widths = [GABA_TEP(s).widths(1:2), GABA_TEP(s).widths(2), GABA_TEP(s).widths(3:5)];
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

% identify individual peaks - baseline datasets
% k = 1; p = 1; m = 1; s = 1; e = 1;
for k = 1:length(GABA_TEP(1).peaks)
    % select seed 
    if numel(seed_electrode) ~= numel(seed_peaks)
        disp('Number of seed electrodes does not correspond to the seed peak distribution!')
    end
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
                [peak_center, data_c, data_c_sub] = track_peak(data, header, time_window, k, GABA_TEP(s), buffer, seed);
                
                % fill in outcome structure
                GABA_TEP_subject(p).latency(m, k) = peak_center;
                
                % append centered data
                statement = ['GABA_' GABA_TEP(s).peaks{k} '_centered(counter, :, :) = data_c;'];
                eval(statement)
                
                % append centered subtracted data
                statement = ['GABA_' GABA_TEP(s).peaks{k} '_subtracted(counter, :, :) = data_c_sub;'];
                eval(statement)
                
                % update row counter
                counter = counter + 1;
            end
        end
    end
end
clear counter a k p m t s e data statement peak_center data_c data_c_sub seed

% average across subjects and sessions, save for letswave
for k = 1:length(GABA_TEP(1).peaks)
    % choose centered data
    statement = ['data_i =  GABA_' GABA_TEP(1).peaks{k} '_centered;'];
    eval(statement)
    
    % average across trials (subject x medication)
    for e = 1:size(data_i, 2)
        for i = 1:size(data_i, 3)
            data(1, e, 1, 1, 1, i) =  mean(squeeze(data_i(:, e, i)));         
            GABA_tracking(k).centered(e, i) = mean(squeeze(data_i(:, e, i)));        
        end
    end
    
    % save for LW
    filename = ['TEP tracked ' target ' '  GABA_TEP(1).peaks{k} ' centered'];
    header.name = filename; 
    header.datasize(6) = size(data, 6);
    span = (1 + buffer) * GABA_TEP(1).widths(k);
    header.xstart = - span/2;
    save([filename '.mat'], 'data');
    save([filename '.lw6'], 'header');    
    clear data
    
    % choose centered subtracted data
    statement = ['data_i =  GABA_' GABA_TEP(1).peaks{k} '_subtracted;'];
    eval(statement)
    
    % average across trials (subject x medication)
    for e = 1:size(data_i, 2)
        for i = 1:size(data_i, 3)
            data(1, e, 1, 1, 1, i) =  mean(squeeze(data_i(:, e, i)));         
            GABA_tracking(k).subtracted(e, i) = mean(squeeze(data_i(:, e, i)));        
        end
    end
    
    % save for LW
    filename = ['TEP tracked ' target ' '  GABA_TEP(1).peaks{k} ' subtracted'];
    header.name = filename; 
    save([filename '.mat'], 'data');
    save([filename '.lw6'], 'header');    
    clear data    
end
clear k e i data_i statement filename span

% choose EOIs based on centered datasets
for k = 1:length(GABA_TEP(1).peaks)
    % choose data, baseline correct
    data =  GABA_tracking(k).centered;
    
    % choose peak polarity
    if strcmp(GABA_TEP(1).peaks{k}(1), 'P')
        operation = 'max';
    elseif strcmp(GABA_TEP(1).peaks{k}(1), 'N')
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
    x_lim = ((1 + buffer) * GABA_TEP(1).widths(k))/2;
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
    title(GABA_TEP(1).peaks{k})
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
function gmfp_topoplot(header, data, x_pos, x_start, map_lims)
varargin = {'maplimits' map_lims 'shading' 'interp' 'whitebk' 'on'};

% fetch data to display
x_visual = ceil((x_pos - x_start)/header.xstep);
vector = data(:, x_visual);

%fetch chanlocs
chanlocs = header.chanlocs;

%parse data and chanlocs 
k=1;
for chanpos=1:size(chanlocs,2);
    vector2(k)=double(vector(chanpos));
    chanlocs2(k)=chanlocs(chanpos);
    k=k+1;
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
end

