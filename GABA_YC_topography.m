%% GABA-AD: TEPs - TOPOGRAPHICAL ANALYSIS
% Written by Dominika for GABA-AD project (2022)
% 
% Colection of scripts to visualize the outcome of the topographical analysis performed in Ragu 
%   --> figures are saved in a folder 'GABA_YC_figures'
% 
% 1) load Ragu output
% 2) load GFP data
% 3) BASELINE: test of topographic consistency
% 4) ALL COMPARISONS: TANOVA
% 5) EFFECT OF ALPRAZOLAM: GFP

%% parameters
clear all; clc

% ----- adjustable parameters -----
% dataset
prefix = 'GABA';
group = 'YC';
participant = 1:20;
medication = {'placebo' 'alprazolam'}; 
time = {'pre' 'post'};
peaks_M1 = {'N17' 'P30' 'N45' 'P60' 'N100' 'P180'};
peaks_AG = {'P25' 'N40' 'P50' 'N75' 'N100' 'P180'};

% visualization (in ms)
time_window = [-50, 300];
analysis_window = [10, 300];
x_delta = 0.5;
shade = 0.2;

% statistics
z = 1.96;
alpha = 0.05;
% --------------------------------

% navigate to the Git folder --> source of saved default files
folder_git = uigetdir(pwd, 'Choose the Git folder');

% visualization 
figure_counter = 1;
x = time_window(1):x_delta:time_window(2);

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

% input/output folders
folder_results = uigetdir(pwd, 'Choose the Results folder');
folder_input = [folder_results '\GABA_' group '_microstates'];
folder_figures = [folder_results '\GABA_' group '_figures'];

%% 1) LOAD RAGU OUTPUT
data = struct;

% baseline - target
load([folder_input '\GABA_' group '_baseline_target.mat'])
data.target = rd;

% baseline - intensity
load([folder_input '\GABA_' group '_baseline_intensity.mat'])
data.intensity = rd;

% medication - M1 TS
load([folder_input '\GABA_' group '_medication_M1-TS.mat'])
data.M1_TS = rd;

% medication - M1 CS
load([folder_input '\GABA_' group '_medication_M1-CS.mat'])
data.M1_CS = rd;

% medication - AG
load([folder_input '\GABA_' group '_medication_AG.mat'])
data.AG = rd;

% set up condiitions
condition = fieldnames(data);
clear rd

%% 2) LOAD GFP
% mean GFP - from all subjects (baseline)
load([folder_results '\GABA_' group '_variables\GABA_' group '_M1_TEPs.mat'], 'GABA_GFP_mean')
GFP_M1 = GABA_GFP_mean;
load([folder_results '\GABA_' group '_variables\GABA_' group '_AG_TEPs.mat'], 'GABA_GFP_mean')
GFP_AG = GABA_GFP_mean;

data_GFP_all(1, :, :, :) = squeeze(GFP_M1(:, :, 2, :));
data_GFP_all(2, :, :, :) = squeeze(GFP_M1(:, :, 1, :));
data_GFP_all(3, :, :, :) = squeeze(GFP_AG(:, :, :, :));

disp(['data_GFP_all - datasize: ' num2str(size(data_GFP_all))])
clear GABA_GFP_mean GFP_AG GFP_M1

% Choose subjects to display - alprazolam condition
for c = 1:3
    statement = ['subj2plot(c, :) = (data.' condition{c + 2} '.IndFeature == 1);'];
    eval(statement)
end

% mean GFP - from selected subjects (alprazolam)
load([folder_results '\GABA_' group '_variables\GABA_' group '_M1_TEPs.mat'], 'GABA_data')
TEP_M1 = GABA_data;
load([folder_results '\GABA_' group '_variables\GABA_' group '_AG_TEPs.mat'], 'GABA_data')
TEP_AG = GABA_data;

data_GFP_selected_subj(1, :, :, :, :) = squeeze(std(TEP_M1(:, :, 2, :, 1:30, :), 1, 5));
data_GFP_selected_subj(2, :, :, :, :) = squeeze(std(TEP_M1(:, :, 1, :, 1:30, :), 1, 5));
data_GFP_selected_subj(3, :, :, :, :) = squeeze(std(TEP_AG(:, :, :, 1:30, :), 1, 4));
disp(['data_GFP_selected_subj - datasize: ' num2str(size(data_GFP_selected_subj))])

data_GFP_selected(1, :, :, :) = squeeze(std(mean(TEP_M1(:, :, 2, subj2plot(1, :), 1:30, :), 4), 1, 5));
data_GFP_selected(2, :, :, :) = squeeze(std(mean(TEP_M1(:, :, 1, subj2plot(2, :), 1:30, :), 4), 1, 5));
data_GFP_selected(3, :, :, :) = squeeze(std(mean(TEP_AG(:, :, subj2plot(3, :), 1:30, :), 3), 1, 4));
disp(['data_GFP_selected - datasize: ' num2str(size(data_GFP_selected))])

clear TEP_M1 TEP_AG GABA_data

%% 3) BASELINE: TOPOGRAPHIC CONSISTENCY
% define time intervals of inconsistency
incons_t = {[], [10:x_delta:21, 33:x_delta:50.5], [10:x_delta:13.5, 26:x_delta:32, 38.5:x_delta:52]};
incons = {zeros(1, length(x)), zeros(1, length(x)), zeros(1, length(x))};
for n = 1:length(incons)
    incons_t{n} = (incons_t{n} - time_window(1))/x_delta;
    incons{n}(incons_t{n}) = 1;
end
clear n incons_t

% plot baseline GFP with TCT results for all stimuli
for c = 1:3
    % prepare baseline data (average sessions)
    data_visual = squeeze(mean(data_GFP_all(c, :, 1, :), 2))';

    % plot GFP + TCT
    fig = plot_TCT(x, data_visual, c, figure_counter, incons)
%     title(['overall GFP: ' condition{c}])

    % name and save figure
    figure_name = ['GFP_incons_' condition{c + 2}];
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')

    % update figure counter
    figure_counter = figure_counter + 1 ;
end
clear c i data_visual fig yl I P lgd figure_name incons 

%% 4) ALL COMPARISONS: TANOVA  
% define comparisons of interest
comps = {1, 1, 3, 3, 3};
for c = 1:length(condition)
    statement = ['comps_name = {data.' condition{c} '.strF1, data.' condition{c} '.strF2, [data.' condition{c} '.strF1 ''_'' data.' condition{c} '.strF2]};'];
    eval(statement)    
    comps{2, c} = comps_name{comps{1, c}};
end
clear c comps_name

% plot TANOVA results for selected comparisons
for c = 1:numel(condition)
    % ----- TANOVA p value -----
    % load data
    statement = ['data_visual = squeeze(data.' condition{c} '.PTanova(1, 1 + comps{1, c}, :, 1));'];
    eval(statement)       

    % plot the p-value timecourse
    fig = plot_p(x, data_visual, figure_counter, analysis_window);                   

    % name and save figure
    figure_name = ['TANOVA_' condition{c} '_' comps{2, c}];
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')

    % update figure counter
    figure_counter = figure_counter + 1 ;       
    
    % ----- explained variance -----
    % load data
    statement = ['data_visual = squeeze(data.' condition{c} '. TExpVar{1}(1, 1 + comps{1, c}, :));']; 
    eval(statement)       

    % plot EV timecourse
    fig = plot_EV(x, data_visual, figure_counter, analysis_window);

    % name and save figure
    figure_name = ['TANOVA_EV_' condition{c} '_' comps{2, c}];
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')

    % update figure counter
    figure_counter = figure_counter + 1 ;       
end
clear c a statement fig figure_name data_visual

% define intervals for t-map extraction
signif_intervals = {};
for c = 1:2
    % define intervals of significance
    statement = ['data_p = squeeze(data.' condition{c} '.PTanova(1, 1 + comps{1, c}, :, 1));'];
    eval(statement)     
    signif = data_p < 0.05;
    
    % define intervals of inconsistency
    switch c
        case 1 
            incons_TANOVA = incons{2} | incons{3};
        case 2
            incons_TANOVA = incons{1} | incons{2};
    end
    
    % extract intervals of interest
    signif(find(incons_TANOVA)) = false;
    lims = [];
    counter = 1;
    for i = 2:length(signif)
        if ~signif(i-1) & signif(i)
            lims(counter, 1) = i * x_delta + time_window(1);
        elseif signif(i-1) & ~signif(i)
            lims(counter, 2) = i * x_delta + time_window(1);
            counter = counter + 1;
        end
    end
    signif_intervals(c) = {lims};
    
end
clear c i statement data_p incons_TANOVA signif lims

%% 5) EFFECT OF ALPRAZOLAM: GFP
n_perms = 1500;
threshold = norminv(1 - alpha/2);
for c = 1:3
    for m = 1:2
        % ----- identify time intervals of significant change in GFP -----
        % create data matrix
        n_subj = size(data_GFP_selected_subj, 4);
        for p = 1:n_subj
            data_diff(:, p) = squeeze(data_GFP_selected_subj(c, m, 2, p, :) - data_GFP_selected_subj(c, m, 1, p, :))';
        end
        data_diff = data_diff(:, subj2plot(c, :));
        n_subj = size(data_diff, 2);
        data_zero = zeros(size(data_diff));
        data_perm = cat(2, data_diff, data_zero);
        clear p data_diff data_zero
        
        % generate true condition labels
        labels = (1:n_subj*2) <= n_subj;
        
        % run the permutation
        diffs_permuted = zeros(length(x), n_perms);
        for p = 1:n_perms
            % shuffle condition label vector
            labels_mixed = labels(randperm(2*n_subj));

            % compute and store difference time series
            mean_1 = mean(data_perm(:, labels_mixed == 1), 2);
            mean_2 = mean(data_perm(:, labels_mixed == 0), 2);
            diffs_permuted(:, p) = mean_1 - mean_2;
        end
        clear p labels_mixed mean_1 mean_2
        
        % compute z-score difference
        diff_observed = mean(data_perm(:, labels == 1), 2);
        z_diff = (diff_observed - mean(diffs_permuted, 2)) ./ std(diffs_permuted, [], 2);
        
        % statistically threshold the final result
        z_thresholded = z_diff;
        z_thresholded(abs(z_thresholded) < threshold) = 0;
        z_GFP_all(1, c, m, :) =  z_thresholded;
        
        % ----- correct for too small clusters -----
        % find cluster sizes under the null hypothesis
        cluster_size = zeros(n_perms, 1);
        for p = 1:n_perms
            % compute z-score difference
            z_diff_perm = (diffs_permuted(:, p) - mean(diffs_permuted, 2)) ./ std(diffs_permuted,[],2);

            % threshold
            z_diff_perm(abs(z_diff_perm) < threshold) = 0;

            % identify clusters
            islands = bwconncomp(logical(z_diff_perm));

            % find cluster sizes
            if length(islands.PixelIdxList) == 0
                cluster_size(p) = 0;
            else
                clustNs = cellfun(@length, islands.PixelIdxList);
                cluster_size(p) = max(clustNs);
            end
        end
        clear p z_diff_perm islands clustNs

        % compute cluster threshold
        cluster_threshold(c, m) = prctile(cluster_size, 100 - alpha*100);

%         % show distribution of cluster sizes
%         figure
%         histogram(cluster_size)
%         xlabel('Cluster size (time points)'), ylabel('Count')

        % find islands
        islands = bwconncomp(logical(z_thresholded));

        % find and remove any subthreshold islands
        z_thresholded_cc = z_thresholded;
        if islands.NumObjects > 0
            for i = 1:islands.NumObjects
                if numel(islands.PixelIdxList{i}) < cluster_threshold(c, m)
                    z_thresholded_cc(islands.PixelIdxList{i}) = 0;
                end
            end
        end
        z_GFP_all(2, c, m, :) =  z_thresholded_cc;
        clear i cluster_size islands
        
        % ----- plot the data -----
        % extract comparison names
        statement = ['comps = {data.' condition{c} '.strF1, data.' condition{c} '.strF2, [data.' condition{c} '.strF1 ''-'' data.' condition{c} '.strF2]};'];
        eval(statement)
    
        % choose the data
        data_visual = squeeze(data_GFP_selected(c, m, :, :));
        
        % plot GFP timecourse
        col = cat(1, [0 0 0], colours(2*(m-1) + 2, :));
        fig = plot_GFP_z(x, data_visual, c, figure_counter, analysis_window, col, z_thresholded_cc);
        
        % name and save figure
        figure_name = ['GFP_z_' condition{c} '_' medication{m}];
        savefig([folder_figures '\' figure_name '.fig'])
        saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')

        % update figure counter
        figure_counter = figure_counter + 1 ;         
    end
end
clear c m p i n_subj data_perm labels diffs_permuted diff_observed z_diff threshold n_perms threshold...
    z_thresholded z_thresholded_cc statement comps data_visual col fig figure_name 

%% 6) GFP - RAGU OUTPUT
% add post-medication intervals of inconsistency
incons_t = {[], [10:x_delta:24, 33.5:x_delta:53], [10:x_delta:63]};
incons = {zeros(1, length(x)), zeros(1, length(x)), zeros(1, length(x))};
for n = 1:length(incons)
    incons_t{n} = (incons_t{n} - time_window(1))/x_delta;
    incons{n}(incons_t{n}) = 1;
end
clear n incons_t

% loop through conditions
for c = 3:length(condition) 
    % plot GFP per medication 
    for m = 1:2
        % extract comparison names
        statement = ['comps = {data.' condition{c} '.strF1, data.' condition{c} '.strF2, [data.' condition{c} '.strF1 ''-'' data.' condition{c} '.strF2]};'];
        eval(statement)
    
        % choose the data
        data_visual = squeeze(data_GFP_selected(c - 2, m, :, :));
        
        % extract intervals of significant effect of time
        comp_n = find(strcmp(comps, 'time')) + 1;
        statement = ['data_p = squeeze(data.' condition{c} '.GFPPTanova(1, comp_n, :, 1));'];
        eval(statement)          
        signif(1, :) = data_p < alpha;
        
        % extract intervals of significant effect of time
        statement = ['data_p = squeeze(data.' condition{c} '.GFPPTanova(1, 4, :, 1));'];
        eval(statement)          
        signif(2, :) = data_p < alpha;
        clear data_p comp_n
        
        % plot GFP timecourse
        col = cat(1, [0 0 0], colours(2*(m-1) + 2, :));
        fig = plot_GFP(x, data_visual, c, figure_counter, analysis_window, col, incons, signif);
        
        % add legend
        lgd = legend({'baseline' ['post ' medication{m}]});
        lgd.FontSize = 18; 
        lgd.Location = 'northoutside';
        lgd.NumColumns = 2;
        
        % name and save figure
        figure_name = ['GFP_' condition{c} '_' medication{m}];
        savefig([folder_figures '\' figure_name '.fig'])
        saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')

        % update figure counter
        figure_counter = figure_counter + 1 ;  
    end
            
    % GFP p value
    for a = 3
        % load data
        statement = ['data_visual = squeeze(data.' condition{c} '.GFPPTanova(1, a+1, :, 1));'];
        eval(statement)       
        
        % plot the p-value timecourse
        fig = plot_p(x, data_visual, figure_counter, analysis_window);        
%         title([condition{c} ' GFP: factor ''' comps{a} ''''])               
                
        % name and save figure
        figure_name = ['GFP_' condition{c} '_' comps{a}];
        savefig([folder_figures '\' figure_name '.fig'])
        saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')

        % update figure counter
        figure_counter = figure_counter + 1 ;       
    end 
    
    % explained variance
    for a = 3
        % load data
        statement = ['data_visual = squeeze(data.' condition{c} '. GFPExpVar{1}(1, a+1, :));']; 
        eval(statement)       
        
        % plot EV timecourse
        fig = plot_EV(x, data_visual, figure_counter, analysis_window);
        
        % name and save figure
        figure_name = ['GFP_EV_' condition{c} '_' comps{a}];
        savefig([folder_figures '\' figure_name '.fig'])
        saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')

        % update figure counter
        figure_counter = figure_counter + 1 ;       
    end  
end
clear c m a statement fig figure_name data_visual col lgd

%% functions
function fig = plot_TCT(x, data_visual, stim, figure_counter, incons)
    % launch the figure
    fig = figure(figure_counter);
    hold on
    
    % set limits of y
    switch stim
        case 1
            yl = [0 4.6];
        case 2
            yl = [0 2.1];
        case 3
            yl = [0 3.1];
    end
    ylim(yl)
    
    % shade interpolated interval 
    rectangle('Position', [-5, yl(1), 15, yl(2) - yl(1)], 'FaceColor', [0.99 0.73 0.73], 'EdgeColor', 'none')

    % plot intervals of inconsistency
    I = area(x, incons{stim} * yl(2));
    I.FaceColor = [0.85 0.85 0.85];
    I.EdgeColor = 'none';
    
    % plot GFP   
    P = plot(x, data_visual, 'Color', [0 0 0], 'LineWidth', 3);
        
    % TMS stimulus
    line([0, 0], yl, 'Color', [0.88 0.08 0.08], 'LineWidth', 3)
    
    % other parameters
    xlabel('time (ms)')
    ylabel('amplitude (\muV)')
    set(gca, 'FontSize', 18) 
    xlim([x(1), x(end)])  
    set(gca, 'layer', 'top');
    hold off
    
    % change figure size
    fig.Position = [500 500 750 300];
end
function fig = plot_p(x, data_visual, figure_counter, analysis_window)
    % launch the figure
    fig = figure(figure_counter);
    hold on
    
    % define intervals of significance
    signif_05 = data_visual < 0.05;
    signif_01 = data_visual < 0.01;
    
    % shade intervals of significance
    I(1) = area(x, signif_05);
    I(1).FaceColor = [1 0.73 0.73];
    I(1).EdgeColor = 'none';
    I(2) = area(x, signif_01);
    I(2).FaceColor = [1 0.44 0.44];
    I(2).EdgeColor = 'none';
    
    % plot p
    P = plot(x, data_visual, 'Color', [0 0 0], 'LineWidth', 3);
    
    % other parameters  
    set(gca, 'xtick', [])
%     xlabel('time (ms)')
    ylabel('probability')
    set(gca, 'FontSize', 18)
    xlim([analysis_window(1) + 1, analysis_window(2)])    
    set(gca, 'layer', 'top');
    
    % change figure size
    fig.Position = [500 500 750 300];
    
    hold off
end
function fig = plot_EV(x, data_visual, figure_counter, analysis_window)
    % launch the figure
    fig = figure(figure_counter);
    hold on
    
    % shade area of explained variance
    A = area(x, data_visual*100);
    A.FaceColor = [0.65 0.65 0.65];
    A.EdgeColor = 'none';

    % other parameters    
%     set(gca, 'xtick', [])
    xlabel('time (ms)')
    ylabel('% explained variance')
    set(gca, 'FontSize', 18)
    xlim([analysis_window(1) + 2, analysis_window(2)])   
    set(gca, 'layer', 'top');

    % change figure size
    fig.Position = [500 500 750 300];
    
    hold off
end
function fig = plot_GFP(x, data_visual, stim, figure_counter, analysis_window, col, incons, signif)
    % launch the figure
    fig = figure(figure_counter);
    hold on

    % set limits of y
    switch stim
        case 3
            yl = [-2 5.5];
        case 4
            yl = [-1.5 2.5];
        case 5
            yl = [-1.5 3.5];
    end
    ylim(yl)
       
    % shade intervals of significance
    for s = 1:2
        % effect of time
        T(s) = area(x, signif(1, :) * yl(s));
        T(s).FaceColor = [0    0.4471    0.7412];
        T(s).FaceAlpha = 0.5;
        T(s).EdgeColor = 'none';
        T(s).Annotation.LegendInformation.IconDisplayStyle = 'off';
        
        % effect of interaction
        N(s) = area(x, signif(2, :) * yl(s));
        N(s).FaceColor = [0.9804    0.5569    0.5569];
        N(s).FaceAlpha = 0.5;
        N(s).EdgeColor = 'none';
        N(s).Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    
    % plot intervals of inconsistency
    for i = 1:2
        I(i) = area(x, incons{stim - 2} * yl(i));
        I(i).FaceColor = [0.85 0.85 0.85];
        I(i).EdgeColor = 'none';
        I(i).Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    
    % plot the GFP course
    for t = 1:size(data_visual, 1)
        plot(x, data_visual(t, :), 'Color', col(t, :), 'LineWidth', 3);
    end
    
    % shade area of GFP difference
    diff = data_visual(2, :) - data_visual(1, :);
    F = fill([x fliplr(x)],[diff zeros(1, length(x))], [0.4 0.4 0.4], 'linestyle', 'none');
    F.Annotation.LegendInformation.IconDisplayStyle = 'off';
%     F(1) = fill([x((diff > 0)) fliplr(x((diff > 0)))],[diff(diff > 0) zeros(1, length(x((diff > 0))))],...
%         [0.8784 0.0784 0.0784], 'linestyle', 'none', 'facealpha', 0.4);
%     F(2) = fill([x((diff < 0)) fliplr(x((diff < 0)))],[diff(diff < 0) zeros(1, length(x((diff < 0))))],...
%         [0 0.4471 0.7412], 'linestyle', 'none', 'facealpha', 0.8);
%     F(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
%     F(2).Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    % add a horizontal line
%     line(x, zeros(1, length(x)), 'color', [0 0 0], 'LineWidth', 2, 'LineStyle', ':')

    % other parameters    
%     set(gca, 'xtick', [])
    xlabel('time (ms)')
    ylabel('amplitude(\muV)')
    set(gca, 'FontSize', 18)
    xlim([analysis_window(1) + 2, analysis_window(2)])   
    set(gca, 'layer', 'top');

    % change figure size
    fig.Position = [500 500 550 400];
    
    hold off
end
function fig = plot_GFP_z(x, data_visual, stim, figure_counter, analysis_window, col, z_thresholded)
    % launch the figure
    fig = figure(figure_counter);
    hold on

    % set limits of y
    switch stim
        case 1
            yl = [-2 5.5];
        case 2
            yl = [-1.5 2.5];
        case 3
            yl = [-1.5 3.5];
    end
    ylim(yl)
           
    % plot intervals of significant difference
    for i = 1:2
        I(i) = area(x, logical(z_thresholded) * yl(i));
        I(i).FaceColor = [0.9882    0.7608    0.7608];
        I(i).EdgeColor = 'none';
        I(i).Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    
    % plot the GFP course
    for t = 1:size(data_visual, 1)
        plot(x, data_visual(t, :), 'Color', col(t, :), 'LineWidth', 3);
    end
    
    % shade area of GFP difference
    diff = data_visual(2, :) - data_visual(1, :);
    F = fill([x fliplr(x)],[diff zeros(1, length(x))], [0.4 0.4 0.4], 'linestyle', 'none');
    F.Annotation.LegendInformation.IconDisplayStyle = 'off';

    % other parameters
    xlabel('time (ms)')
    ylabel('amplitude(\muV)')
    set(gca, 'FontSize', 18)
    xlim([analysis_window(1) + 2, analysis_window(2)])   
    set(gca, 'layer', 'top');

    % change figure size
    fig.Position = [500 500 550 400];
    
    hold off
end
