%% GABA-AD: TEPs VISUALIZATION
% Written by Dominika for GABA-AD project (2020-21)
% 
% Colection of scripts to visualy explore extracted TMS-EEG variables



%% parameters
clear all; clc

% ----- adjustable parameters -----
% dataset
prefix = 'GABA';
group = 'YC';
participant = 1:20;
medication = {'placebo' 'alprazolam'}; 
time = {'pre' 'post'};
stimulus = {'M1 - CS' 'M1 - TS' 'AG' 'M1 - ppTMS'};
peaks_M1 = {'N17' 'P30' 'N45' 'P60' 'N100' 'P180'};
peaks_AG = {'P25' 'N40' 'P50' 'N75' 'N100' 'P180'};

% visualization (in ms)
time_window = [-50, 300];
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

% load default header
load([folder_git '\GABA_header_default.mat'])
labels = {header.chanlocs.labels};

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

% output folders
folder_results = uigetdir(pwd, 'Choose the Results folder');
folder_figures = [folder_results '\GABA_' group '_figures'];

%% 1) LOAD DATA
% extracted outcome variables
load([folder_results '\GABA_' group '_statistics\GABA_YC_results.mat'])

% individual TEPs
load([folder_results '\GABA_' group '_variables\GABA_' group '_M1_TEPs.mat'], 'GABA_data')
TEP_M1 = GABA_data;
disp(['M1 TEPs - datasize: ' num2str(size(TEP_M1))])

load([folder_results '\GABA_' group '_variables\GABA_' group '_AG_TEPs.mat'], 'GABA_data')
TEP_AG = GABA_data;
disp(['AG TEPs - datasize: ' num2str(size(TEP_AG))])

TEP.data(1, :, :, :, :, :) = squeeze(TEP_M1(:, :, 1, :, 1:30, :));
TEP.data(2, :, :, :, :, :) = squeeze(TEP_M1(:, :, 2, :, 1:30, :));
TEP.data(3, :, :, :, :, :) = squeeze(TEP_AG(:, :, :, 1:30, :));
TEP.data(4, :, :, :, :, :) = squeeze(TEP_M1(:, :, 3, :, 1:30, :));
clear GABA_data_mean TEP_AG TEP_M1

% calculate mean values
for s = 1:size(TEP.data, 1)
    for m = 1:length(medication)
        for t = 1:length(time)
            for e = 1:size(TEP.data, 5)
                for i = 1:size(TEP.data, 6)
                    TEP.mean(s, m, t, e, i) = mean(squeeze(TEP.data(s, m, t, :, e, i)));
                    TEP.SD(s, m, t, e, i) = std(squeeze(TEP.data(s, m, t, :, e, i)));
                    TEP.SEM(s, m, t, e, i) = TEP.SD(s, m, t, e, i)/sqrt(length(participant));
                    TEP.CI(s, m, t, e, i) = TEP.SEM(s, m, t, e, i) * z;
                end
            end  
        end 
    end
end
clear s m t e i 

% SICI
load([folder_results '\GABA_' group '_variables\GABA_' group '_M1_TEPs.mat'], 'GABA_SICI')

%% ) BASELINE TEPs
% plot the baseline butterfly plots for all stimuli
for s = 1:3
    % prepare baseline data data (placebo session)
    data_visual = squeeze(TEP.mean(s, 1, 1, :, :));

    % launch the figure
    fig = figure(figure_counter);
    hold on

    % set limits of y
    switch s
        case 1
            yl = [-4.2 4.2];
        case 2
            yl = [-9 9];
        case 3
            yl = [-5 5];
    end
    ylim(yl)

    % shade interpolated interval 
    rectangle('Position', [-5, yl(1) + 0.02, 15, yl(2) - yl(1) - 0.02], 'FaceColor', [0.99 0.73 0.73], 'EdgeColor', 'none')

    % loop through channels to plot
    for c = 1:size(data_visual, 1)     
        P(c) = plot(x, data_visual(c, :), 'Color', [0.65 0.65 0.65], 'LineWidth', 1.5);
    end
    
    % Cz electrode
    P(end+1) =  plot(x, data_visual(18, :), 'Color', [0 0 0], 'LineWidth', 3);
    
    % TMS stimulus
    line([0, 0], yl, 'Color', [0.88 0.08 0.08], 'LineWidth', 3)

    % other parameters
    title(['baseline TEP: ' stimulus{s}])
    xlabel('time (ms)')
    ylabel('amplitude (\muV)')
    set(gca, 'FontSize', 14)
    xlim([x(1), x(end)])    
    lgd = legend(P(end), 'Cz electrode', 'Box', 'off');
    lgd.FontSize = 14;
    lgd.Position = [0.225 0.05 0.4 0.3];
    hold off

    % name and save figure
    if length(stimulus{s}) > 2
        figure_name = ['TEP_' stimulus{s}(1:2) '_' stimulus{s}(end-1:end) '_baseline'];
    else
        figure_name = ['TEP_' stimulus{s}(1:2) '_baseline'];
    end
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')

    % update figure counter
    figure_counter = figure_counter + 1 ;
end
clear s data_visual fig yl c P lgd figure_name

%% ) EFFECT OF ALPRAZOLAM: M1 TEPs
% ----- decide output parameters -----
electrode = {'C3'};
y_limits = [-2.5 4.5; -8 8; -3.5 6];
% ------------------------------------
% identify electrode(s) to plot
for e = 1:length(electrode)
    electrode_n(e) = find(strcmp(labels, electrode{e}));
end
clear e

% loop through stimuli and medication
for s = 1:3   
    for m = 1:length(medication)
        % prepare data
        data_visual = squeeze(TEP.mean(s, m, :, electrode_n, :));
        SEM_visual = squeeze(TEP.SEM(s, m, :, electrode_n, :));
        
        % choose colours
        col(1, :) = [0 0 0];
        col(2, :) = colours((m-1)*2 + 2, :);
        
        % plot the figure
        [fig, figure_counter] = plot_TEP(data_visual, SEM_visual, x, y_limits(s, :), col, figure_counter);
        
        % name and save figure
        if length(stimulus{s}) > 2
            figure_name = sprintf('TEP_%s_%s_%s_%s', stimulus{s}(1:2), stimulus{s}(end-1:end), medication{m}, electrode{1});
        else
            figure_name = sprintf('TEP_%s_%s_%s', stimulus{s}(1:2), medication{m}, electrode{1});
        end
        savefig([folder_figures '\' figure_name '.fig'])
        saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')        
    end
end
clear s m data_visual SEM_visual col

% plot butterfly plots
for s = 1:3   
    for m = 1:length(medication)
        % prepare data
        data_visual = squeeze(TEP.mean(s, m, :, :, :));
        
        % choose colours
        col(1, :) = [0.65 0.65 0.65];
        col(2, :) = colours((m-1)*2 + 2, :);
        
        % plot the figure
        [fig, figure_counter] = plot_TEP_butterfly(data_visual, x, y_limits(s, :), col, figure_counter);
        
        % name and save figure
        if length(stimulus{s}) > 2
            figure_name = sprintf('TEP_%s_%s_%s_butterfly', stimulus{s}(1:2), stimulus{s}(end-1:end), medication{m});
        else
            figure_name = sprintf('TEP_%s_%s_butterfly', stimulus{s}(1:2), medication{m});
        end
        savefig([folder_figures '\' figure_name '.fig'])
        saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')        
    end
end
clear s m data_visual SEM_visual col

%% ) EFFECT OF ALPRAZOLAM: M1 TEP PEAKS
% calculate group mean values
for m = 1:length(medication)
    for s = 1:2
        for k = 1:length(peaks_M1)
            avg_amp(m, s, k) = mean(GABA_YC_results.TEP_M1(m).amplitude.change(:, s, k));
            avg_amp_sem(m, s, k) = std(GABA_YC_results.TEP_M1(m).amplitude.change(:, s, k))/sqrt(length(participant));
        end
    end
end
clear m s k 
disp(['Datasize: ' num2str(size(avg_amp))])

% plot a barplot 
for s = 1:2
    % choose peaks to include
    if s == 1
        peak_n = 2:6;
    else
        peak_n = 1:6;
    end     

    % get data
    data_visual = squeeze(avg_amp(:, s, peak_n))';
    sem_visual = squeeze(avg_amp_sem(:, s, peak_n))';
    
    % launch the figure
    fig = figure(figure_counter);
    hold on
    barplot = bar(data_visual, 'EdgeColor', 'none');
    for b = 1:size(data_visual, 2)
        barplot(b).FaceColor = colours((b-1)*2 + 2, :);
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
%     title(sprintf('M1, %s: change in TEP amplitude', stimulus{s}))
    ylabel('\Delta amplitude (\muV \pmSEM)');
    xlabel('TEP component')
    set(gca, 'xtick', 1:length(peak_n), 'xticklabel', peaks_M1(peak_n))
    set(gca, 'Fontsize', 14)
    legend(barplot, medication, 'Location', 'southwest', 'fontsize', 14)

    % name and save figure
    figure_name = sprintf('TEP_%s_change_barplot', stimulus{s});
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.svg'])
    figure_counter = figure_counter + 1;
end
clear s peak_n m k fig barplot b ngroups nbars groupwidth i x_bar figure_name avg_amp avg_amp_sem

% plot a boxplot
for s = [1 2]
    % determine peaks to plot
    if s == 1
        peaks_visual = peaks_M1(2:end);
    else 
        peaks_visual = peaks_M1;
    end
    
    % get individual data 
    data_amplitude = [];
    if s == 1 
        for m = 1:length(medication)       
            for k = 1:length(peaks_visual)
                data_amplitude(m, k, :) = GABA_YC_results.TEP_M1(m).amplitude.change(:, s, k+1);
            end
        end
    else
        for m = 1:length(medication)       
            for k = 1:length(peaks_visual)
                data_amplitude(m, k, :) = GABA_YC_results.TEP_M1(m).amplitude.change(:, s, k);
            end
        end
    end
    clear m k

    % plot 
    label = ['TEP peak' peaks_visual];
    fig = plot_box(data_amplitude, 'amplitude', medication, colours([2 4], :), label, figure_counter)
    hold off

    % save figure 
    fig_name = sprintf('TEP_%s_change_boxplot', stimulus{s});
    savefig([folder_figures '\' fig_name '.fig'])
    saveas(fig, [folder_figures '\' fig_name  '.svg'])  
    figure_counter = figure_counter + 1;
end

%% ) EFFECT OF ALPRAZOLAM: TEP - MEP CORRELATION
% get the data
for m = 1:length(medication)
    for p = 1:length(participant)    
        % TEP change = x 
        x_corr((m-1)*length(participant) + p, :) = squeeze(GABA_YC_results.TEP_M1(m).amplitude.change(p, 2, :));
        
        % MEP change = y
        y_corr((m-1)*length(participant) + p) = GABA_YC_results.MEP(m).amplitude.change(p, 1);
        
        % choose colours
        marker_col((m-1)*length(participant) + p, :) = colours((m-1)*2 + 2, :);
    end
end
clear m p

% calculate and plot for each component of interest
stats_corr = table;
stats_coefs = struct;
for p = 1:length(peaks_M1)
    % choose x and y
    data_corr(:, 1) = x_corr(:, p);
    data_corr(:, 2) = y_corr';
    
    % preliminary plot
    figure(figure_counter)
    scatter(data_corr(:, 1), data_corr(:, 2))
    xlabel(peaks_M1{p})
    ylabel('MEP')
    figure_counter = figure_counter + 1;
    
    % ----- compute the correlation ----        
    % prepare linear model: y ~ 1 + x
    data_model = fitlm(data_corr(:, 1), data_corr(:, 2), 'VarNames', {peaks_M1{p} 'MEP'});
    [cor_coef, cor_p] = corr(data_corr, 'Type', 'Pearson');
    stats_coefs.r(p) = cor_coef(1, 2);
    stats_coefs.p(p) = cor_p(1, 2);
    clear cor_coef cor_p
    
    % extract statistics
    data_table = data_model.Coefficients;
    data_table.R2(1) = data_model.Rsquared.Ordinary;
    data_table.R2(2) = data_model.Rsquared.Adjusted;
    data_table.Properties.RowNames = {};
    data_table.Variable(1) = {'(intercept)'};
    data_table.Variable(2) = data_model.VariableNames(1);
    data_table = data_table(:, [end, 1:end - 1]);
    stats_corr = [stats_corr; data_table];
    
    % ------------- plot -------------
    % plot data + regression line
    fig = figure(figure_counter);
    hold on
    plot_corr(data_model, data_corr, marker_col, 'Pearson')
    
    % save the figure   
%     figure_name = ['corr_TEPxMEP' peaks_M1{p}];
%     figure_name = 'corr_N17xMEP_change';
%     savefig([folder_figures '\' figure_name '.fig'])
%     saveas(fig, [folder_figures '\' figure_name '.svg'])
    
    % update figure counter
    figure_counter = figure_counter + 1;
end
clear p x_corr y_corr marker_col data_corr data_corr_ranked data_model_ranked data_table fig figure_name

%% ) EFFECT OF ALPRAZOLAM: TEP - RS-EEG CORRELATION
% get the data
data_x = []; data_y = [];
for m = 1:length(medication)
    for p = 1:length(participant)    
        % TEP change = x 
        data_x((m-1)*length(participant) + p, :) = squeeze(GABA_YC_results.TEP_M1(m).amplitude.change(p, 2, 1));
        
        % MEP change = y
        data_y(1, (m-1)*length(participant) + p) = GABA_YC_results.rsEEG(m).sigma.change(p);
        data_y(2, (m-1)*length(participant) + p) = GABA_YC_results.rsEEG(m).delta.change(p);
        data_y(3, (m-1)*length(participant) + p) = GABA_YC_results.rsEEG(m).AAC.change(p);
        data_y(4, (m-1)*length(participant) + p) = GABA_YC_results.rsEEG(m).SE.change.closed(p, 1);
         
        % choose colours
        marker_col((m-1)*length(participant) + p, :) = colours((m-1)*2 + 2, :);
    end
end
data_y = data_y';
clear m p

% get variable names
rsEEG_vars = fieldnames(GABA_YC_results.rsEEG)';
rsEEG_vars = rsEEG_vars(end - size(data_y, 2) + 1 : end);

% calculate and plot for each component of interest
stats_corr = table;
stats_corr_ranked = table;
corr_coef = struct;
for c = 1:size(data_y, 2)
    % choose x and y
    data_corr(:, 1) = data_x;
    data_corr(:, 2) = data_y(:, c);
    
    % preliminary plot
    figure(figure_counter)
    scatter(data_corr(:, 1), data_corr(:, 2))
    xlabel('N17')
    ylabel(rsEEG_vars{c})
    figure_counter = figure_counter + 1;
    
    % ----- compute linear correlation ----     
    % prepare linear model: y ~ 1 + x
    data_model = fitlm(data_corr(:, 1), data_corr(:, 2), 'VarNames', {'N17' rsEEG_vars{c}});
    
    % identify r and p    
    [corr_coef_i, corr_p_i] = corr(data_corr, 'Type', 'Pearson');
    corr_coef.linear.r(c) = corr_coef_i(1, 2);
    corr_coef.linear.p(c) = corr_p_i(1, 2);
    
    % extract statistics
    data_table = data_model.Coefficients;
    data_table.R2(1) = data_model.Rsquared.Ordinary;
    data_table.R2(2) = data_model.Rsquared.Adjusted;
    data_table.Properties.RowNames = {};
    data_table.Variable(1) = {'(intercept)'};
    data_table.Variable(2) = data_model.VariableNames(1);
    data_table = data_table(:, [end, 1:end - 1]);
    stats_corr = [stats_corr; data_table];
    
    % plot data + regression line
    fig = figure(figure_counter);
    hold on
    plot_corr(data_model, data_corr, marker_col, 'Pearson')

    % save the figure   
%     figure_name = ['corr_N17x' rsEEG_vars{c} '_linear'];
%     savefig([folder_figures '\' figure_name '.fig'])
%     saveas(fig, [folder_figures '\' figure_name '.svg'])
    
    % update figure counter
    figure_counter = figure_counter + 1;
       
    % ----- compute ranked correlation ----  
    % rank the data
    for a = 1:size(data_corr, 2)
        [temp, data_corr_ranked(:, a)]  = ismember(data_corr(:, a), unique(data_corr(:, a)));
    end
    clear a temp
    
    % prepare linear model: y ~ 1 + x
    data_model_ranked = fitlm(data_corr_ranked(:, 1), data_corr_ranked(:, 2), 'VarNames', {'N17' rsEEG_vars{c}});
    
    % identify r and p    
    [corr_coef_i, corr_p_i] = corr(data_corr_ranked, 'Type', 'Spearman');
    corr_coef.ranked.r(c) = corr_coef_i(1, 2);
    corr_coef.ranked.p(c) = corr_p_i(1, 2);
    
    % extract statistics
    data_table_ranked = data_model_ranked.Coefficients;
    data_table_ranked.R2(1) = data_model_ranked.Rsquared.Ordinary;
    data_table_ranked.R2(2) = data_model_ranked.Rsquared.Adjusted;
    data_table_ranked.Properties.RowNames = {};
    data_table_ranked.Variable(1) = {'(intercept)'};
    data_table_ranked.Variable(2) = data_model_ranked.VariableNames(1);
    data_table_ranked = data_table_ranked(:, [end, 1:end - 1]);
    stats_corr_ranked = [stats_corr_ranked; data_table_ranked];
    
    % plot data + regression line
    fig = figure(figure_counter);
    hold on
    plot_corr(data_model_ranked, data_corr_ranked, marker_col, 'Spearman')

    % save the figure   
%     figure_name = ['corr_N17x' rsEEG_vars{c} '_ranked'];
%     savefig([folder_figures '\' figure_name '.fig'])
%     saveas(fig, [folder_figures '\' figure_name '.svg'])
    
    % update figure counter
    figure_counter = figure_counter + 1;
end
clear c data_corr data_corr_ranked data_model data_model_ranked...
    data_table data_table_ranked fig figure_name

% plot separately per medication
stats_corr_partial = table;
corr_coef_partial = struct;
for c = 1:size(data_y, 2)
    % choose x and y
    data_corr(:, 1) = data_x;
    data_corr(:, 2) = data_y(:, c);

    % prepare linear models
    data_model_placebo = fitlm(data_corr(1:20, 1), data_corr(1:20, 2), 'VarNames', {'N17' rsEEG_vars{c}});
    data_model_alprazolam = fitlm(data_corr(21:40, 1), data_corr(21:40, 2), 'VarNames', {'N17' rsEEG_vars{c}});
    
    % get coefficients
    [corr_r, corr_p] = corr(data_corr(1:20, :), 'type', 'Pearson');
    corr_coef_partial.placebo.r(c) = corr_r(1, 2); corr_coef_partial.placebo.p(c) = corr_p(1, 2);   
    [corr_r, corr_p] = corr(data_corr(21:40, :), 'type', 'Pearson');
    corr_coef_partial.alprazolam.r(c) = corr_r(1, 2); corr_coef_partial.alprazolam.p(c) = corr_p(1, 2); 
    clear corr_r corr_p

    % plot data + regression line
    fig = figure(figure_counter);
    hold on
    plot_corr_partial(data_model_placebo, data_corr(1:20, :), marker_col(1:20, :), 'Pearson')
    plot_corr_partial(data_model_alprazolam, data_corr(21:40, :), marker_col(21:40, :), 'Pearson')

    % add annotations
    text_pos = [0.90 0.78];
    T(1) = text(0.55, text_pos(1), sprintf( 'y = %1.3f * x', data_model_alprazolam.Coefficients.Estimate(2)),...
        'Units', 'Normalized', 'Color', marker_col(21, :));
    T(2) = text(0.55, text_pos(2), sprintf('y = %1.3f * x', data_model_placebo.Coefficients.Estimate(2)),...
        'Units', 'Normalized', 'Color', marker_col(1, :));
    set(T, 'fontsize', 20, 'fontangle', 'italic', 'fontweight', 'bold'); 

    % save the figure   
%     figure_name = ['corr_N17x' rsEEG_vars{c} '_partial'];
%     savefig([folder_figures '\' figure_name '.fig'])
%     saveas(fig, [folder_figures '\' figure_name '.svg'])

    % update figure counter
    figure_counter = figure_counter + 1;

    % extract statistics
    for m = 1:length(medication)
        statement = ['data_table = data_model_' medication{m} '.Coefficients;'];
        eval(statement)
        statement = ['data_table.R2(1) = data_model_' medication{m} '.Rsquared.Ordinary;'];
        eval(statement)
        statement = ['data_table.R2(2) = data_model_' medication{m} '.Rsquared.Adjusted;'];
        eval(statement)
        data_table.Properties.RowNames = {};
        data_table.Variable(1) = {'(intercept)'};
        data_table.Variable(2) = data_model_placebo.VariableNames(2);
        data_table = data_table(:, [end, 1:end - 1]);
        stats_corr_partial = [stats_corr_partial; data_table];
    end
end
clear x y rsEEG_vars marker_col data_corr data_model_placebo data_model_alprazolam...
    data_table fig text_pos T figure_name m statement

%% ) EFFECT OF ALPRAZOLAM: N17 - BETA CORRELATION
% get the data
for m = 1:length(medication)
    for p = 1:length(participant)
        for t = 1:length(time)
            % determine row
            row = (m-1)*length(participant)*2 + (t-1)*length(participant) + p;
            
            % TEP change = x 
            statement = ['x(row, :) = squeeze(GABA_YC_results.TEP_M1(m).amplitude.' time{t} '(p, 2, 1));'];
            eval(statement)

            % MEP change = y
            statement = ['y(row, :) = squeeze(GABA_YC_results.rsEEG(m).sigma.' time{t} '(p));'];
            eval(statement)

            % choose colours
            marker_col(row, :) = [0.75 0.75 0.75];
        end
    end
end
clear m p t row statement 
    
% ----- compute the correlation ----  
% rank the data
[temp, x_ranked]  = ismember(x, unique(x));
[temp, y_ranked]  = ismember(y, unique(y));
clear a temp

% prepare linear model: y ~ 1 + x
data_model = fitlm(x_ranked, y_ranked, 'VarNames', {'N17' 'low beta'});
    
% extract statistics
data_table = data_model.Coefficients;
data_table.R2(1) = data_model.Rsquared.Ordinary;
data_table.R2(2) = data_model.Rsquared.Adjusted;
data_table.Properties.RowNames = {};
data_table.Variable(1) = {'(intercept)'};
data_table.Variable(2) = data_model.VariableNames(1);

% ------------- plot -------------
% plot data + regression line
fig = figure(figure_counter);
hold on
plot_corr(data_model, cat(2, x_ranked, y_ranked), marker_col, 'Spearman')
    
% save the figure   
figure_name = 'corr_N17xbeta_all';
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.svg'])
    
% update figure counter
figure_counter = figure_counter + 1;
clear p x_corr y_corr marker_col data_corr data_corr_ranked data_model_ranked data_table fig figure_name

%% ) EFFECT OF ALPRAZOLAM: AG TEPs
peaks_SICI = {'N17' 'P60' 'N100'}; 
% calculate group mean values
for m = 1:length(medication)
    for k = 1:length(peaks_AG)
        avg_amp(m, k) = mean(GABA_YC_results.TEP_AG(m).amplitude.change(:, k));
        avg_amp_sem(m, k) = std(GABA_YC_results.TEP_AG(m).amplitude.change(:, k))/sqrt(length(participant));
    end
end
clear m k 
disp(['Datasize: ' num2str(size(avg_amp))]) 

% get data
data_visual = avg_amp';
sem_visual = avg_amp_sem';
    
% launch the figure
fig = figure(figure_counter);
hold on
barplot = bar(data_visual, 'EdgeColor', 'none');
for b = 1:size(data_visual, 2)
    barplot(b).FaceColor = colours((b-1)*2 + 2, :);
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
% title('AG: change in TEP amplitude')
ylabel('\Delta amplitude (\muV \pmSEM)');
xlabel('TEP component')
set(gca, 'xtick', 1:length(peaks_AG), 'xticklabel', peaks_AG)
set(gca, 'Fontsize', 14)
legend(barplot, medication, 'Location', 'southwest', 'fontsize', 14)

% name and save figure
figure_name = 'TEP_AG_amplitude_change';
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.svg'])

% update figure counter
figure_counter = figure_counter + 1;
clear fig barplot b ngroups nbars groupwidth i x_bar figure_name avg_amp avg_amp_sem

% plot a boxplot
for s = 3
    % determine peaks to plot
    peaks_visual = peaks_AG;
    
    % get individual data 
    data_amplitude = [];
    for m = 1:length(medication)       
        for k = 1:length(peaks_visual)
            data_amplitude(m, k, :) = GABA_YC_results.TEP_AG(m).amplitude.change(:, k);
        end
    end
    clear m k

    % plot 
    label = ['TEP peak' peaks_visual];
    fig = plot_box(data_amplitude, 'amplitude', medication, colours([2 4], :), label, figure_counter)
    hold off

    % save figure 
    fig_name = sprintf('TEP_%s_change_boxplot', stimulus{s});
    savefig([folder_figures '\' fig_name '.fig'])
    saveas(fig, [folder_figures '\' fig_name  '.svg'])  
    figure_counter = figure_counter + 1;
end

%% ) BASELINE SICI: AMPLITUDE CHANGE
% ----- decide output parameters ----- 
outliers = [2];
% ------------------------------------
% identify participants without outliers
WO_idx = find(~ismember(participant, outliers));

% get data 
answer = questdlg('What do you wish to plot?', 'Choose data', 'baseline', 'both', 'baseline'); 
switch answer
    case 'baseline'
        % subset data
        for m = 1:length(medication) 
            for k = 1:length(peaks_M1) 
                data_visual(k, m) = mean(GABA_YC_results.SICI(m).TEP.pre(WO_idx, k));  
                sem_visual(k, m) = std(GABA_YC_results.SICI(m).TEP.pre(WO_idx, k))/sqrt(length(WO_idx)); 
            end 
        end 
        % choose colours
        col = colours([2, 4], :); 
        
    case 'both'
        % subset data
        for m = 1:length(medication)
            for t = 1:length(time)
                for k = 1:length(peaks_M1) 
                    statement = ['data_visual(k, (m-1)*2 + t) = mean(GABA_YC_results.SICI(m).TEP.' time{t} '(WO_idx, k));'];  
                    eval(statement)
                    statement = ['sem_visual(k, (m-1)*2 + t) = std(GABA_YC_results.SICI(m).TEP.' time{t} '(WO_idx, k))/sqrt(length(WO_idx));']; 
                    eval(statement)
                end 
            end
        end 
        % choose colours
        col = colours([2, 1, 4, 3], :); 
end
clear answer t m k statement
         
% ----- plot a barplot -----
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
title('SICI: change in TEP peak amplitude') 
ylabel('\Delta amplitude (\muV \pmSEM)'); 
xlabel('TEP component') 
set(gca, 'xtick', 1:6, 'xticklabel', peaks_M1) 
set(gca, 'Fontsize', 14) 
legend(barplot, medication, 'Location', 'southeast', 'fontsize', 14) 
legend off
 
% name and save figure 
figure_name = 'TEP-SICI_change_both'; 
savefig([folder_figures '\' figure_name '.fig']) 
saveas(fig, [folder_figures '\' figure_name '.svg']) 
 
% update figure counter, clear 
figure_counter = figure_counter + 1; 
clear m k data_visual sem_visual fig barplot col b ngroups nbars groupwidth i x_bar 

% ----- plot a boxplot -----
% determine peaks to plot
peaks_visual = peaks_M1;
    
% get individual data 
for m = 1:length(medication)       
    for k = 1:length(peaks_visual)
        data_amplitude(m, k, :) = GABA_YC_results.SICI(m).TEP.pre(:, k);
    end
end
        
% plot 
label = ['TEP peak' peaks_visual];
fig = plot_box(data_amplitude, 'amplitude', medication, colours([2 4], :), label, figure_counter)
hold off

% save figure 
fig_name = 'TEP-SICI_change_boxplot';
savefig([folder_figures '\' fig_name '.fig'])
saveas(fig, [folder_figures '\' fig_name  '.svg'])  
figure_counter = figure_counter + 1;

%% ) BASELINE SICI: TEP - MEP CORRELATION
% ----- decide output parameters -----
peaks_SICI = {'N17' 'P60' 'N100'}; 
outliers = [2];
% ------------------------------------
% identify participants without outliers
WO_idx = find(~ismember(participant, outliers));

% get the data
x_corr = [], y_corr = [];
for m = 1:length(medication)
    for p = 1:length(WO_idx)
        % TEP SICI = x 
        x_corr(:, (m-1)*length(WO_idx) + p) = GABA_YC_results.SICI(m).TEP.pre(WO_idx(p), contains(peaks_M1, peaks_SICI))';
        
        % MEP SICI = y
        y_corr((m-1)*length(WO_idx) + p) = 100 - GABA_YC_results.SICI(m).MEP.pre(WO_idx(p));
        
        % choose colours
        marker_col((m-1)*length(WO_idx) + p, :) = colours((m-1)*2 + 2, :);
    end
end
clear m p

% calculate and plot for each component of interest
stats_corr = table;
for p = 1:length(peaks_SICI)
    % ----- compute the correlation -----
    % choose x and y
    data_corr(:, 1) = x_corr(p, :)'*-1;
    data_corr(:, 2) = y_corr';
    
    % rank the data
    for a = 1:size(data_corr, 2)
        [temp, data_corr_ranked(:, a)]  = ismember(data_corr(:, a), unique(data_corr(:, a)));
    end
    clear a temp
    
    % prepare linear model: y ~ 1 + x
    data_model_ranked = fitlm(data_corr_ranked(:, 1), data_corr_ranked(:, 2), 'VarNames', {[peaks_SICI{p} ' SICI'] 'MEP SICI'});
    
    % extract statistics
    data_table = data_model_ranked.Coefficients;
    data_table.R2(1) = data_model_ranked.Rsquared.Ordinary;
    data_table.R2(2) = data_model_ranked.Rsquared.Adjusted;
    data_table.Properties.RowNames = {};
    data_table.Variable(1) = {'(intercept)'};
    data_table.Variable(2) = data_model_ranked.VariableNames(1);
    data_table = data_table(:, [end, 1:end - 1]);
    stats_corr = [stats_corr; data_table];
    
    % ------------- plot -------------
    % plot data + regression line
    fig = figure(figure_counter);
    hold on
    plot_corr(data_model_ranked, data_corr_ranked, marker_col, 'Spearman')

    % save the figure   
%     figure_name = ['corr_SICI-MEPxSICI-TEP_' peaks_SICI{p}];
%     savefig([folder_figures '\' figure_name '.fig'])
%     saveas(fig, [folder_figures '\' figure_name '.svg'])
    
    % update figure counter
    figure_counter = figure_counter + 1;
end
clear p x_corr y_corr marker_col data_corr data_corr_ranked data_model_ranked data_table fig figure_name

% for p = 1:length(peaks_SICI)
%     % fit
%     [curvefit,gof,output] = fit(x, y, 'exp1');
%     
%     % plot
%     figure(figure_counter)
%     subplot(2, 2, 1)
%     plot(curvefit,x,y)
%     
%     subplot(2, 2, 2)
%     scatter(x,log(y))  
%     
%     subplot(2, 2, 3)
%     plot(curvefit,x,y,'Residuals')
%     
%     subplot(2, 2, 4)
%     plot(curvefit,x,y,'predfunc')
%     
%     % compute coefficients
%     n = length(x);
%     y2 = log(y);
%     j = sum(x);
%     k = sum(y2);
%     r2 = sum(x .* y2);
%     m = sum(y2.^2);
%     c = f.b * (r2 - j * k / n);
%     d = m - k^2 / n;
%     corr = sqrt(c/d);
%     std_err = sqrt((d - c) / (n - 2)); 
%     
% end
% clear x_corr y_corr marker_col p x y curvefit n y2 r2 j k m c d corr std_err gof output

%% ) BASELINE SICI: SUBTRACTION
% ----- decide output parameters -----
electrode = {'C3'}; 
outliers = [2];
% ------------------------------------
% identify participants without outliers
WO_idx = find(~ismember(participant, outliers));

% load default labels
load([folder_git '\GABA_header_default.mat'])
labels = {header.chanlocs.labels};
clear header

% time axis
x = [-50:0.5:300];

% calculate mean data, save for letswave
data = [];
counter = 1;
for m = 1:2
    for p = 1:20
        for e = 1:30
            for i = 1:701
                data(counter, e, 1, 1, 1, i) = GABA_SICI.individual(m, 1, p, e, i);
            end
        end
        counter = counter + 1;
    end
end
for e = 1:30
    for i = 1:701
        data_mean(e, i) = mean(GABA_SICI.individual(:, 1, :, e, i), [1, 3]);
        data_sem(e, i) = std(squeeze(data(:, e, 1, 1, 1, i)), 0, 1)/sqrt(length(participant));
    end
end
clear m p e i counter
save('merged SICI baseline.mat', 'data')

% save header
header = lwdata.header;
header.name = 'merged SICI baseline';
header.datasize = size(data);
header.xstart = -0.05;
header.chanlocs = lwdata.header.chanlocs(1:30);
save('merged SICI baseline.lw6', 'header')

% plot SICI for each chosen electrode
for e = 1:length(electrode)
%     % prepare baseline data data at C3 (placebo session)
%     e_n = find(strcmp(labels, electrode{e}));
%     data_visual = squeeze(GABA_SICI.mean(1, 1, e_n, :))';
%     CI_visual = squeeze(GABA_SICI.SEM(1, 1, e_n, :))';

    % prepare baseline data data at C3 
    e_n = find(strcmp(labels, electrode{e}));
    data_visual = squeeze(data_mean(e_n, :))';
    CI_visual = squeeze(data_sem(e_n, :))';

    % launch the figure
    fig = figure(figure_counter);
    hold on

    % set limits of the figure
    plot(x, data_visual + CI_visual, 'b:', 'LineWidth', 0.5)
    plot(x, data_visual - CI_visual, 'b:', 'LineWidth', 0.5)
    yl = get(gca, 'ylim'); yl(1) = yl(1) - 0.2; yl(2) = yl(2) + 0.3;
    cla, hold on
    
    % shade interpolated interval 
    rectangle('Position', [-5, yl(1) + 0.02, 15, yl(2) - yl(1) - 0.02], 'FaceColor', [0.99 0.73 0.73], 'EdgeColor', 'none')

    % plot data       
    P = plot(x, data_visual, 'Color', [0 0 0], 'LineWidth', 3);
    F = fill([x fliplr(x)],[data_visual + CI_visual fliplr(data_visual - CI_visual)], ...
        [0 0 0], 'FaceAlpha', 0.25, 'linestyle', 'none');

    % mark TMS stimulus and zero line
    line([0, 0], yl, 'Color', [0.88 0.08 0.08], 'LineWidth', 2)
    line(x, zeros(1, length(x)), 'Color', [0.88 0.08 0.08], 'LineWidth', 2, 'LineStyle', ':')
    
    % add other parameters
    xlabel('time (ms)')
    ylabel('amplitude (\muV \pm SEM)')
    set(gca, 'FontSize', 14)
    xlim([x(1), x(end)])
    ylim(yl)
    
    % add legend
    lgd = legend(P, [electrode{e} ' electrode'], 'Box', 'off');
    lgd.FontSize = 14;
    lgd.Position = [0.55 0.05 0.4 0.3];
    hold off
    
    % save figure
    figure_name = ['TEP_SICI_' electrode{e} '_baseline'];
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.svg'])
end

% update figure counter
figure_counter = figure_counter + 1 ;

clear electrode e e_n x data_visual CI_visual fig yl P F lgd figure_name 

%% ) BASELINE SICI vs. EFFECT OF ALPRAZOLAM 
% ----- decide output parameters -----
peaks_SICI = {'N17' 'P60' 'N100'}; 
% ------------------------------------
stats_corr = [];
corr_coef = struct;
for a = 1:length(peaks_SICI)
    % get the data - only from the alprazolam session
    x_corr = []; y_corr = [];
    for p = 1:length(participant)    
        % SICI of N17 = x 
        x_corr(p) = GABA_YC_results.SICI(2).TEP.pre(p, contains(peaks_M1, peaks_SICI{a}))';

        % alprazolam induced change in N17 = y
        y_corr(p) = GABA_YC_results.TEP_M1(2).amplitude.change(p, 2, contains(peaks_M1, peaks_SICI{a}));

        % choose colours
        marker_col(p, :) = colours(4, :);
    end
    clear p
    if figure_counter == 1
        figure_counter = 2
    end
    figure(1)
    scatter(x_corr, y_corr)
    
    % prepare the data & linear model based on the type of correlation
    data_corr = [];
    data_corr(:, 1) = x_corr';
    data_corr(:, 2) = y_corr';
    corr_type = questdlg('What type of correlation do you want to run?', 'Correlatio', 'Pearson', 'Spearman', 'Pearson'); 
    
    % recover  r and p
    [corr_r, corr_p] = corr(data_corr, 'type', corr_type);
    corr_coef.r(a) = corr_r(1, 2); corr_coef.p(a) = corr_p(1, 2);  
    clear corr_r corr_p
    
    % rank data if necessary
    if strcmp(corr_type, 'Spearman')
        [temp, data_corr(:, 1)]  = ismember(data_corr(:, 1), unique(data_corr(:, 1)));
        [temp, data_corr(:, 2)]  = ismember(data_corr(:, 2), unique(data_corr(:, 2)));
        clear temp
    end
    
    % prepare the linear model
    data_model = fitlm(data_corr(:, 1), data_corr(:, 2),...
        'VarNames', {sprintf('%s SICI', peaks_SICI{a}) sprintf('%s  modulation', peaks_SICI{a})});
    
    % extract statistics
    data_table = data_model.Coefficients;
    data_table.R2(1) = data_model.Rsquared.Ordinary;
    data_table.R2(2) = data_model.Rsquared.Adjusted;
    data_table.Properties.RowNames = {};
    data_table.Variable(1) = {'(intercept)'};
    data_table.Variable(2) = data_model.VariableNames(1);
    data_table = data_table(:, [end, 1:end - 1]);
    stats_corr = [stats_corr; data_table];

    % plot data + regression line
    fig = figure(figure_counter);
    hold on
    plot_corr(data_model, data_corr, marker_col, corr_type)

    % save the figure   
%     figure_name = sprintf('corr_SICI_alprazolam_%s_%s', peaks_SICI{a}, corr_type);
%     savefig([folder_figures '\' figure_name '.fig'])
%     saveas(fig, [folder_figures '\' figure_name '.svg'])

    % update figure counter
    figure_counter = figure_counter + 1;
end
clear a x_corr y_corr marker_col data_corr data_model fig figure_name

%% ) BASELINE SICI vs. EFFECT OF ALPRAZOLAM: N17 vs. RS-EEG
% ----- decide output parameters -----
rsEEG_vars = {'sigma' 'delta' 'AAC' 'SE'}; 
% ------------------------------------
% get the data - only from the alprazolam session
for p = 1:length(participant)    
    % TEP SICI = x 
    x_corr(p) = GABA_YC_results.SICI(2).TEP.pre(p, contains(peaks_M1, 'N17'))';

    % RS-EEG = y
    for v = 1:length(rsEEG_vars)
        if v < 4
            statement = ['y_corr(v, p) = GABA_YC_results.rsEEG(2).' rsEEG_vars{v} '.change(p);'];
            eval(statement)
        else
            y_corr(v, p) = GABA_YC_results.rsEEG(2).SE.change.closed(p);
        end
    end

    % choose colours
    marker_col(p, :) = colours(4, :);
end
clear p v statement

% calculate and plot for each component of interest
stats_corr = table;
corr_coef = struct;
for v = 1:length(rsEEG_vars)
    % ----- compute the correlation -----
    % choose x and y
    figure(figure_counter)
    data_corr = [];
    data_corr(:, 1) = x_corr';
    data_corr(:, 2) = y_corr(v, :)';
    scatter(data_corr(:, 1), data_corr(:, 2))
    
    % update figure counter
    figure_counter = figure_counter + 1;
    
    % recover  r and p
    [corr_r, corr_p] = corr(data_corr, 'type', 'Pearson');
    corr_coef.r(v) = corr_r(1, 2); corr_coef.p(v) = corr_p(1, 2);  
    clear corr_r corr_p

    % prepare linear model: y ~ 1 + x
    data_model = fitlm(data_corr(:, 1), data_corr(:, 2), 'VarNames', {'N17 SICI' rsEEG_vars{v}});
    
    % extract statistics
    data_table = data_model.Coefficients;
    data_table.R2(1) = data_model.Rsquared.Ordinary;
    data_table.R2(2) = data_model.Rsquared.Adjusted;
    data_table.Properties.RowNames = {};
    data_table.Variable(1) = {'(intercept)'};
    data_table.Variable(2) = data_model.VariableNames(1);
    data_table = data_table(:, [end, 1:end - 1]);
    stats_corr = [stats_corr; data_table];
    
    % ------------- plot -------------
    % plot data + regression line
    fig = figure(figure_counter);
    hold on
    plot_corr(data_model, data_corr, marker_col, 'Pearson')

    % save the figure   
%     figure_name = ['corr_SICI-N17_' rsEEG_vars{v}];
%     savefig([folder_figures '\' figure_name '.fig'])
%     saveas(fig, [folder_figures '\' figure_name '.svg'])

    % update figure counter
    figure_counter = figure_counter + 1;
end
clear v x_corr y_corr marker_col data_corr data_model fig figure_name data_table

%% ) MEP CHANGE
% extract data
data_visual = cat(2, GABA_YC_results.MEP(1).amplitude.change(:, 1), GABA_YC_results.MEP(2).amplitude.change(:, 1));
label = ['medication' medication];

% plot the boxplot
fig = plot_scatter(data_visual, figure_counter, colours([2, 4], :), 'comparison', label)

% save figure 
fig_name = 'MEP_change_boxplot';
savefig([folder_figures '\' fig_name '.fig'])
saveas(fig, [folder_figures '\' fig_name  '.svg'])  
figure_counter = figure_counter + 1;

%% ) MEP-SICI CHANGE
% extract data
data_visual = cat(2, GABA_YC_results.SICI(1).MEP.post', GABA_YC_results.SICI(2).MEP.post');
label = ['medication' medication];

% plot the boxplot
fig = plot_scatter(data_visual, figure_counter, colours([2, 4], :), 'comparison', label)

% save figure 
fig_name = 'SICI_MEP_post_boxplot';
savefig([folder_figures '\' fig_name '.fig'])
saveas(fig, [folder_figures '\' fig_name  '.svg'])  
figure_counter = figure_counter + 1;

%% FUNCTIONS
function [fig, figure_counter] = plot_TEP(data_visual, SEM_visual, x, y_limits, col, figure_counter)
    % launch the figure
    fig = figure(figure_counter);
    hold on

    % shade interpolated interval 
    plot(x, data_visual, 'b:', 'LineWidth', 0.5);
    rectangle('Position', [-5, y_limits(1), 13, y_limits(2) - y_limits(1)], 'FaceColor', [0.99 0.73 0.73], 'EdgeColor', 'none')

    % loop through datasets to plot
    for t = 1:size(data_visual, 1)        
        P(t) = plot(x, data_visual(t, :), 'Color', col(t, :), 'LineWidth', 2.5);
        F(t) = fill([x fliplr(x)],[data_visual(t, :) + SEM_visual(t, :) fliplr(data_visual(t, :) - SEM_visual(t, :))], ...
            col(t, :), 'FaceAlpha', 0.2, 'linestyle', 'none');
    end

    % TMS stimulus
    line([0, 0], y_limits, 'Color', [0.88 0.08 0.08], 'LineWidth', 3, 'LineStyle', '--')

    % other parameters
    xlabel('time (ms)')
    ylabel('amplitude (\muV)')
    set(gca, 'FontSize', 16)
    ylim(y_limits)
    xlim([x(1), x(end)])  
    set(gca, 'Layer', 'Top')
    
    % update figure counter
    figure_counter = figure_counter + 1;
end
function [fig, figure_counter] = plot_TEP_butterfly(data_visual, x, y_limits, col, figure_counter)
    % launch the figure
    fig = figure(figure_counter);
    hold on

    % shade interpolated interval 
    rectangle('Position', [-5, y_limits(1), 13, y_limits(2) - y_limits(1)], 'FaceColor', [0.99 0.73 0.73], 'EdgeColor', 'none')

    % loop through datasets to plot
    for t = 1:size(data_visual, 1)   
        data = squeeze(data_visual(t, :, :));
        for e = 1:size(data, 1)
            P((t-1)*size(data, 1) + e) = plot(x, data(e, :), 'Color', col(t, :), 'LineWidth', 1.5);
        end
    end

    % TMS stimulus
    line([0, 0], y_limits, 'Color', [0.88 0.08 0.08], 'LineWidth', 3, 'LineStyle', '--')

    % other parameters
    xlabel('time (ms)')
    ylabel('amplitude (\muV)')
    set(gca, 'FontSize', 16)
    ylim(y_limits)
    xlim([x(1), x(end)])  
    set(gca, 'Layer', 'Top')
    
    % update figure counter
    figure_counter = figure_counter + 1;
end
function plot_corr(data_model, data_corr, marker_col, corr_type)
% plot correlation
plot_cor = plotAdded(data_model);

% identify p    
[cor_coef, cor_p] = corr(data_corr, 'Type', corr_type);

% adjust parameters    
set(gca, 'FontSize', 20)
xlabel(data_model.VariableNames{1}); 
ylabel(data_model.VariableNames{2});
plot_cor(2).Color = [0 0 0]; plot_cor(2).LineWidth = 4; 
plot_cor(3).Color = [0 0 0]; plot_cor(3).LineWidth = 2;
legend off
title('')

% replot markers
for c = 1:size(data_corr, 1)
    scatter(data_corr(c, 1), data_corr(c, 2), 90, marker_col(c, :), 'filled');
    hold on
end

% add annotations
text_pos = [0.90 0.78 0.66];
% rectangle('Position', [2, 27, 18, 12], 'FaceColor',[1 1 1], 'EdgeColor', [0.5 0.5 0.5])
T(1) = text(0.1, text_pos(1), sprintf( 'y = %1.3f * x', data_model.Coefficients.Estimate(2)), 'Units', 'Normalized');
T(2) = text(0.1, text_pos(2), sprintf('R^2 = %1.3f', data_model.Rsquared.Adjusted), 'Units', 'Normalized');
T(3) = text(0.1, text_pos(3), sprintf('p = %1.5f', cor_p(1, 2)), 'Units', 'Normalized');
set(T(1), 'fontsize', 20, 'fontangle', 'italic'); 
set(T(2), 'fontsize', 20, 'fontweight', 'bold'); 
set(T(3), 'fontsize', 20, 'fontweight', 'bold'); 
end
function plot_corr_partial(data_model, data_corr, marker_col, corr_type)
% plot correlation
plot_cor = plotAdded(data_model);

% adjust parameters    
set(gca, 'FontSize', 20)
xlabel(data_model.VariableNames{1}); 
ylabel(data_model.VariableNames{2});
plot_cor(2).Color = [0 0 0]; plot_cor(2).LineWidth = 4; 
plot_cor(3).Color = [0 0 0]; plot_cor(3).LineWidth = 2;
legend off
title('')

% replot markers
for c = 1:size(data_corr, 1)
    scatter(data_corr(c, 1), data_corr(c, 2), 90, marker_col(c, :), 'filled');
    hold on
end
end
function fig = plot_scatter(data_visual, figure_counter, col, varargin)
    % determine compared variables
    b = find(strcmpi(varargin, 'comparison'));
    if ~isempty(b)
        comparison = varargin{b + 1};
    end
    
    % launch the figure
    fig = figure(figure_counter);
    hold on
    
    % plot the markers
    data_group = cat(1, ones(size(data_visual, 1), 1), 2 * ones(size(data_visual, 1), 1));
    data_y = cat(1, data_visual(:, 1), data_visual(:, 2));
    beeswarm(data_group, data_y)
    
    % add boxplot
    for a = 1:2
        boxchart(a * ones(size(data_visual, 1), 1), data_visual(:, a), 'BoxFaceColor', col(a, :))
    end
    
    % reorder layers
    child = get(gca, 'Children');
    set(gca,'Children',[child(3) child(4) child(1) child(2)])
    
    % font
    set(gca, 'Fontsize', 18)
    
    % x label
    xlabel(comparison{1})
    set(gca, 'xtick', 1:size(data_visual, 2), 'xticklabel', comparison(2:end))
    
    % y label
    ylabel('\Delta amplitude (% baseline)')
    ylim([0, 120])
    
    % add zero line
    xl = get(gca, 'xlim');
    line(xl, [100, 100], 'LineStyle', ':', 'LineWidth', 2, 'Color', [0 0 0])
end
function beeswarm(x,y,varargin)
% function xbee = beeswarm(x,y)
%
% Input arguments:
%   x               column vector of groups (only tested for integer)
%   y               column vector of data
%
% Optional input arguments:
%   sort_style      ('nosort' - default | 'up' | 'down' | 'fan' | 'rand' | 'square' | 'hex')
%   corral_style    ('none' default | 'gutter' | 'omit' | 'rand')
%   dot_size        relative. default=1
%   overlay_style   (false default | 'box' | 'sd' | 'ci')
%   use_current_axes (false default | true)
%   colormap        (lines default | 'jet' | 'parula' | 'r' | Nx3 matrix of RGB values]
%
% Output arguments:
%   xbee            optimized layout positions
%
% Known Issues:
%       x locations depend on figure aspect ratio. resizing the figure window and rerunning may give different results
%       setting corral to 'none' still has a gutter when the width is large
%
% Usage example:
% 	x = round(rand(150,1)*5);
%   y = randn(150,1);
%   beeswarm(x,y,3,'sort_style','up','overlay_style','ci')
%
% % Ian Stevenson, CC-BY 2019

p = inputParser;
addRequired(p,'x')
addRequired(p,'y')
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addOptional(p,'sort_style','nosort')
addOptional(p,'corral_style','none')
addOptional(p,'dot_size',20/sqrt(length(x)),validScalarPosNum)
addOptional(p,'overlay_style',false)
addOptional(p,'use_current_axes',false)
addOptional(p,'colormap','lines')
addOptional(p,'MarkerFaceColor','')
addOptional(p,'MarkerFaceAlpha',1)
addOptional(p,'MarkerEdgeColor', 'black')
parse(p,x,y,varargin{:});

% extra parameters
rwid = .05; % width of overlay box/dash

dcut=8; % spacing factor
nxloc=512; % resolution for optimization
chanwid = .9; % percent width of channel to use
yl = [min(y) max(y)]; % default y-limits
asp_rat = 1;
keep_hold = false;

% get aspect ratio for a figure window
if isfinite(p.Results.dot_size)
    if ~p.Results.use_current_axes
        % make new axes
        s=scatter(x,y);
        xl=[min(x)-.5 max(x)+.5];
    else
        xl=xlim();
    end
    yl=ylim();
    pasp_rat = get(gca,'PlotBoxAspectRatio');
    dasp_rat = get(gca,'DataAspectRatio');
    asp_rat = pasp_rat(1)/pasp_rat(2);
    
    % pix-scale
    pf = get(gcf,'Position');
    pa = get(gca,'Position');
    as = pf(3:4).*pa(3:4); % width and height of panel in pixels
    dcut = dcut*sqrt(p.Results.dot_size)/as(1)*(range(unique(x))+1);
    cla
end

% sort/round y for different plot styles
yorig=y;
switch lower(p.Results.sort_style)
    case 'up'
        [y,sid]=sort(y);
    case 'fan'
        [~,sid]=sort(abs(y-mean(y)));
        sid=[sid(1:2:end); sid(2:2:end)];
        y=y(sid);
    case 'down'
        [y,sid]=sort(y,'descend');
    case 'rand'
        sid=randperm(length(y));
        y=y(sid);
    case 'square'
        nxloc=.9/dcut;
%         [~,e,b]=histcounts(y,ceil((range(x)+1)*chanwid*nxloc/2/asp_rat));
        edges = linspace(min(yl),max(yl),ceil((range(x)+1)*chanwid*nxloc/asp_rat));
        [~,e,b]=histcounts(y,edges);
        y=e(b)'+mean(diff(e))/2;
        [y,sid]=sort(y);
    case 'hex'
        nxloc=.9/dcut;
%         [~,e,b]=histcounts(y,ceil((range(x)+1)*chanwid*nxloc/2/sqrt(1-.5.^2)/asp_rat));
        edges = linspace(min(yl),max(yl),ceil((range(x)+1)*chanwid*nxloc/sqrt(1-.5.^2)/asp_rat));
        [n,e,b]=histcounts(y,edges);
        oddmaj=0;
        if sum(mod(n(1:2:end),2)==1)>sum(mod(n(2:2:end),2)==1),
            oddmaj=1;
        end
        y=e(b)'+mean(diff(e))/2;
        [y,sid]=sort(y);
        b=b(sid);
    otherwise
        sid=1:length(y);
end
x=x(sid);
yorig=yorig(sid);
[ux,~,ic] = unique(x);
% rmult=(range(ux)+1)*2;
rmult=5;

% for each group...
for i=1:length(ux)
    fid = find(ic==i);   
    
    % set of possible x locations
    xi = linspace(-chanwid/2*rmult,chanwid/2*rmult,nxloc*rmult+(mod(nxloc*rmult,2)==0))'+ux(i);

    % rescale y to that things are square visually
    zy=(y(fid)-min(yl))/(max(yl)-min(yl))/asp_rat*(range(ux)+1)*chanwid;
    
    % precalculate y distances so that we only worry about nearby points
    D0=squareform(pdist(zy))<dcut*2;    
    
    if length(fid)>1
        % for each data point in the group sequentially...
        for j=1:length(fid)
            if strcmp(lower(p.Results.sort_style),'hex')
                xi = linspace(-chanwid/2*rmult,chanwid/2*rmult,nxloc*rmult+(mod(nxloc*rmult,2)==0))'+ux(i);
                if mod(b(fid(j)),2)==oddmaj
                    xi = linspace(-chanwid/2*rmult,chanwid/2*rmult,nxloc*rmult+(mod(nxloc*rmult,2)==0))'+ux(i)+mean(diff(xi))/2;
                end
            end
            zid = D0(j,1:j-1);
            e = (xi-ux(i)).^2; % cost function
            if ~strcmp(lower(p.Results.sort_style),'hex') && ~strcmp(lower(p.Results.sort_style),'square')
                if sum(zid)>0
                    D = pdist2([xi ones(length(xi),1)*zy(j)], [x(fid(zid)) zy(zid)]);
                    D(D<=dcut)=Inf;
                    D(D>dcut & isfinite(D))=0;
                    e = e + sum(D,2) + randn(1)*10e-6; % noise to tie-break
                end
            else
                if sum(zid)>0
                    D = pdist2([xi ones(length(xi),1)*zy(j)], [x(fid(zid)) zy(zid)]);
                    D(D==0)=Inf;
                    D(D>dcut & isfinite(D))=0;
                    e = e + sum(D,2) + randn(1)*10e-6; % noise to tie-break
                end
            end

            if strcmp(lower(p.Results.sort_style),'one')
                e(xi<ux(i))=Inf;
            end
            [~,mini] = min(e);
            if mini==1 && rand(1)>.5, mini=length(xi); end
            x(fid(j)) = xi(mini);
        end
    end
%     x(fid)=x(fid)-median(x(fid))+ux(i); % center x locations by median
end

if strcmp(lower(p.Results.sort_style),'randn')
    x=ux(ic)+randn(size(ic))/4;
end

% corral any points outside of the channel
out_of_range = abs(x-ux(ic))>chanwid/2;
switch lower(p.Results.corral_style)
    case 'gutter'
        id = (x-ux(ic))>chanwid/2;
        x(id)=chanwid/2+ux(ic(id));
        id = (x-ux(ic))<-chanwid/2;
        x(id)=-chanwid/2+ux(ic(id));
    case 'omit'
        x(out_of_range)=NaN;
    case 'random'
        x(out_of_range)=ux(ic(out_of_range))+rand(sum(out_of_range),1)*chanwid-chanwid/2;
end

% plot groups and add overlay
if isfinite(p.Results.dot_size)
    if isnumeric(p.Results.colormap)
        cmap=p.Results.colormap;
    else
        cmap = feval(p.Results.colormap,length(ux));
    end
    for i=1:length(ux)
        if isempty(p.Results.MarkerFaceColor')
            scatter(x(ic==i),y(ic==i),p.Results.dot_size*36,'filled','MarkerFaceAlpha',p.Results.MarkerFaceAlpha,'MarkerEdgeColor',p.Results.MarkerEdgeColor,'MarkerFaceColor',cmap(i,:))
        else
            scatter(x(ic==i),y(ic==i),p.Results.dot_size*36,'filled','MarkerFaceAlpha',p.Results.MarkerFaceAlpha,'MarkerEdgeColor',p.Results.MarkerEdgeColor,'MarkerFaceColor',p.Results.MarkerFaceColor)
        end
        hold on
        iqr = prctile(yorig(ic==i),[25 75]);
        switch lower(p.Results.overlay_style)
            case 'box'
                rectangle('Position',[ux(i)-rwid iqr(1) 2*rwid iqr(2)-iqr(1)],'EdgeColor','k','LineWidth',2)
                line([ux(i)-rwid ux(i)+rwid],[1 1]*median(yorig(ic==i)),'LineWidth',3,'Color',cmap(i,:))
            case 'sd'
                line([1 1]*ux(i),mean(yorig(ic==i))+[-1 1]*std(yorig(ic==i)),'Color',cmap(i,:),'LineWidth',2)
                line([ux(i)-2*rwid ux(i)+2*rwid],[1 1]*mean(yorig(ic==i)),'LineWidth',3,'Color',cmap(i,:))
            case 'ci'
                line([1 1]*ux(i),mean(yorig(ic==i))+[-1 1]*std(yorig(ic==i))/sqrt(sum(ic==i))*tinv(0.975,sum(ic==i)-1),'Color',cmap(i,:),'LineWidth',2)
                line([ux(i)-2*rwid ux(i)+2*rwid],[1 1]*mean(yorig(ic==i)),'LineWidth',3,'Color',cmap(i,:))
        end
        
    end
    hold on
    xlim(xl)
    ylim(yl)
end

% unsort so that output matches the original y data
x(sid)=x;
end
function fig = plot_box(data, datatype, condition, col, label, figure_counter)
    % launch the figure
    fig = figure(figure_counter);
    hold on

    % determine x limits
    xl = [0.25, size(data, 2)+0.5];

    % add zero line
    line(xl, [0, 0], 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 1)

    % plot data per condition
    for c = 1:length(condition)  
        % subset data
        data_visual = squeeze(data(c, :, :))';

        % determine x positions
        if c == 1
            data_x(c, :) = (1:size(data_visual, 2))-0.17;
        else
            data_x(c, :) = (1:size(data_visual, 2))+0.17;
        end

        % boxplot
        for t = 1:size(data_visual, 2)
            P(c, t) = boxchart(data_x(c, t) * ones(size(data_visual, 1), 1), data_visual(:, t), ...
                'BoxFaceColor', col(c, :), 'BoxWidth', 0.3, 'WhiskerLineColor', col(c, :), 'MarkerColor', col(c, :));
        end       
    end

    % add legend
    lgd = legend(P(:, 1), condition, 'Location', 'southwest');
    lgd.FontSize = 14;
    legend('boxoff')

    % y label
    if strcmp(datatype, 'amplitude')
        ylabel('\Delta amplitude (\muV \pm SEM)')
    elseif strcmp(datatype, 'latency')
        ylabel('\Delta latency (ms \pm SEM')
    end

    % other parameters
    xlim(xl)
    xlabel(label{1})
    set(gca, 'xtick', 1:size(data_visual, 2), 'xticklabel', label(2:end))
    set(gca, 'FontSize', 14) 
    set(gca, 'layer', 'top');
end