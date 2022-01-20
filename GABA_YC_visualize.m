%% GABA-AD: TMS-EVOKED POTENTIALS 
% Written by Dominika for GABA-AD project (2020-21)
% 
% Colection of scripts to visualy explore extracted TMS-EEG variables: 



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

% visualization 
time_window = [-0.05, 0.3];
shade = 0.2;

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
folder_figures = [folder_results '\GABA_' group '_figures'];

%% 1) LOAD DATA
% extracted outcome variables
load([folder_results '\GABA_' group '_statistics\GABA_YC_results.mat'])

% TEPs
load([folder_results '\GABA_' group '_variables\GABA_' group '_M1_TEPs.mat'], 'GABA_data_mean')
TEP_M1 = GABA_data_mean;
disp(['M1 TEPs - datasize: ' num2str(size(TEP_M1))])

load([folder_results '\GABA_' group '_variables\GABA_' group '_AG_TEPs.mat'], 'GABA_data_mean')
TEP_AG = GABA_data_mean;
disp(['AG TEPs - datasize: ' num2str(size(TEP_AG))])

TEP_data(1, :, :, :, :) = squeeze(TEP_M1(:, :, 1, 1:30, :));
TEP_data(2, :, :, :, :) = squeeze(TEP_M1(:, :, 2, 1:30, :));
TEP_data(3, :, :, :, :) = squeeze(TEP_AG(:, :, 1:30, :));
TEP_data(4, :, :, :, :) = squeeze(TEP_M1(:, :, 3, 1:30, :));

clear GABA_data_mean TEP_AG TEP_M1

%% ) BASELINE TEPs
% time axis
x = [-50:0.5:300];

% plot the baseline butterfly plots for all stimuli
for s = 1:3
    % prepare baseline data data (placebo session)
    data_visual = squeeze(TEP_data(s, 1, 1, :, :));

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
    saveas(fig, [folder_figures '\' figure_name '.png'])

    % update figure counter
    figure_counter = figure_counter + 1 ;
end
clear s data_visual fig yl c P lgd figure_name

%% ) EFFECT OF ALPRAZOLAM - M1 TEPs
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

% plot a barplot for all peaks together
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
    title(sprintf('M1, %s: change in TEP amplitude', stimulus{s}))
    ylabel('\Delta amplitude (\muV \pmSEM)');
    xlabel('TEP component')
    set(gca, 'xtick', 1:length(peak_n), 'xticklabel', peaks_M1(peak_n))
    set(gca, 'Fontsize', 14)
    legend(barplot, medication, 'Location', 'southwest', 'fontsize', 14)

    % name and save figure
    figure_name = ['TEP_M1_' stimulus{s} '_amplitude_change'];
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.png'])

    % update figure counter
    figure_counter = figure_counter + 1;
end
clear s peak_n m k fig barplot b ngroups nbars groupwidth i x_bar figure_name avg_amp avg_amp_sem

%% ) EFFECT OF ALPRAZOLAM - AG TEPs
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
title('AG: change in TEP amplitude')
ylabel('\Delta amplitude (\muV \pmSEM)');
xlabel('TEP component')
set(gca, 'xtick', 1:length(peaks_AG), 'xticklabel', peaks_AG)
set(gca, 'Fontsize', 14)
legend(barplot, medication, 'Location', 'southwest', 'fontsize', 14)

% name and save figure
figure_name = 'TEP_AG_amplitude_change';
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.png'])

% update figure counter
figure_counter = figure_counter + 1;
clear fig barplot b ngroups nbars groupwidth i x_bar figure_name avg_amp avg_amp_sem

%% ) BASELINE SICI
% get data
for m = 1:length(medication)
    for k = 1:length(peaks_M1)
        data_visual(k, m) = mean(GABA_YC_results.SICI(m).TEP.pre(:, k)); 
        sem_visual(k, m) = std(GABA_YC_results.SICI(m).TEP.pre(:, k))/sqrt(length(participant));
    end
end
        
% launch the figure
fig = figure(figure_counter);
hold on
barplot = bar(data_visual, 'EdgeColor', 'none');
col = colours([2, 4], :);
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

% name and save figure
figure_name = 'TEP_SICI_amplitude_baseline';
savefig([folder_figures '\' figure_name '.fig'])
saveas(fig, [folder_figures '\' figure_name '.png'])

% update figure counter, clear
figure_counter = figure_counter + 1;
clear m k data_visual sem_visual fig barplot col b ngroups nbars groupwidth i x_bar

%% ) CORRELATION: TEP SICI x MEP SICI 
% ----- decide output parameters -----
peaks_SICI = {'N17' 'P60' 'N100'};
% ------------------------------------
% variable names
varnames = [peaks_SICI, {'MEP-SICI'}];

% extract TEP data
for m = 1:length(medication)
    for p = 1:length(participant)               
        data_cor((m-1)*length(participant) + p, :) = GABA_YC_results.SICI(m).TEP.pre(p, contains(peaks_M1, peaks_SICI));
    end
end

% extract MEP data
for m = 1:length(medication)
    for p = 1:length(participant)               
        data_cor((m-1)*length(participant) + p, length(peaks_SICI) + 1) = GABA_YC_results.SICI(m).MEP.pre(p);
    end
end

% Bonferroni correction of alpha
alpha_cor = alpha/length(peaks_SICI);

% ----- significant linear correlation -----
% identify significant cases        
[cor_coef, cor_p] = corrcoef(data_cor);
[row, col] = find(cor_p < alpha_cor);

% plot significant cases
for a = 1:length(row)    
    % prepare linear model: y ~ 1 + x
    data_model = fitlm(data_cor(:, row(a)), data_cor(:, col(a)), 'VarNames', [varnames(row(a)) varnames(col(a))]);

    % choose only correlations that show TEP-MEP interactions
    if ismember(row(a), 1:length(peaks_SICI)) && col(a) == length(peaks_SICI) + 1            
        % plot data + regression line
        fig = figure(figure_counter);
        hold on
        plot_cor = plotAdded(data_model);

        % adjust parameters    
        title(sprintf('Linear correlation - baseline:\n%s ~ %s', varnames{col(a)}, varnames{row(a)}))
        xlabel(['change in ' varnames{row(a)}]); ylabel(varnames{col(a)});
        set(gca, 'FontSize', 14)
        plot_cor(1).Marker = 'o'; plot_cor(1).MarkerSize = 8; 
        plot_cor(1).MarkerEdgeColor = colours(2, :); plot_cor(1).MarkerFaceColor = colours(2, :);
        plot_cor(2).Color = colours(3, :); plot_cor(2).LineWidth = 2; 
        plot_cor(3).Color = colours(3, :); plot_cor(3).LineWidth = 2;
        legend off
        if data_model.Coefficients.Estimate(2) > 0
            text_pos = [0.95 0.85 0.75];
        else
            text_pos = [0.25 0.15 0.05];
        end
        T(1) = text(0.05, text_pos(1), sprintf( 'y = %1.3f * x', data_model.Coefficients.Estimate(2)), 'Units', 'Normalized');
        T(2) = text(0.05, text_pos(2), sprintf('R^2 = %1.3f', data_model.Rsquared.Ordinary), 'Units', 'Normalized');
        T(3) = text(0.05, text_pos(3), sprintf('r = %1.3f, p = %1.5f', cor_coef(row(a), col(a)), cor_p(row(a), col(a))), 'Units', 'Normalized');
        set(T(1), 'fontsize', 14, 'fontweight', 'bold', 'fontangle', 'italic'); 
        set(T(2), 'fontsize', 14); 
        set(T(3), 'fontsize', 14, 'color', colours(3, :)); 

        % save figure and continue
        fig_name = ['corr_TEP-SICIxMEP-SICI_linear_' varnames{row(a)} '_' varnames{col(a)}];
        savefig([folder_figures '\' fig_name '.fig'])
        saveas(fig, [folder_figures '\' fig_name  '.png'])  
        figure_counter = figure_counter + 1;
    end
end
% ----- significant ranked correlation -----  

% identify significant cases        
[cor_coef, cor_p] = corr(data_cor, 'Type', 'Spearman');
[row, col] = find(cor_p < alpha_cor);

% rank the data
for a = 1:size(data_cor, 2)
    [temp, data_cor(:, a)]  = ismember(data_cor(:, a), unique(data_cor(:, a)));
end

% plot significant cases
for a = 1:length(row)    
    % prepare linear model: y ~ 1 + x
    data_model = fitlm(data_cor(:, row(a)), data_cor(:, col(a)), 'VarNames', [varnames(row(a)) varnames(col(a))]);

    % choose only correlations that show TEP-MEP interactions
    if ismember(row(a), 1:length(peaks_SICI)) && col(a) == length(peaks_SICI) + 1              
        % plot data + regression line
        fig = figure(figure_counter);
        hold on
        plot_cor = plotAdded(data_model);

        % adjust parameters    
        title(sprintf('Spearman rank correlation - baseline:\n%s ~ %s', varnames{col(a)}, varnames{row(a)}))
        xlabel(['change in ' varnames{row(a)}]); ylabel(varnames{col(a)});
        set(gca, 'FontSize', 14)
        plot_cor(1).Marker = 'o'; plot_cor(1).MarkerSize = 8; 
        plot_cor(1).MarkerEdgeColor = colours(2, :); plot_cor(1).MarkerFaceColor = colours(2, :);
        plot_cor(2).Color = colours(3, :); plot_cor(2).LineWidth = 2; 
        plot_cor(3).Color = colours(3, :); plot_cor(3).LineWidth = 2;
        legend off
        if data_model.Coefficients.Estimate(2) > 0
            text_pos = [0.95 0.85 0.75];
        else
            text_pos = [0.25 0.15 0.05];
        end
        T(1) = text(0.05, text_pos(1), sprintf( 'y = %1.3f * x', data_model.Coefficients.Estimate(2)), 'Units', 'Normalized');
        T(2) = text(0.05, text_pos(2), sprintf('R^2 = %1.3f', data_model.Rsquared.Ordinary), 'Units', 'Normalized');
        T(3) = text(0.05, text_pos(3), sprintf('r = %1.3f, p = %1.5f', cor_coef(row(a), col(a)), cor_p(row(a), col(a))), 'Units', 'Normalized');
        set(T(1), 'fontsize', 14, 'fontweight', 'bold', 'fontangle', 'italic'); 
        set(T(2), 'fontsize', 14); 
        set(T(3), 'fontsize', 14, 'color', colours(3, :)); 

        % save figure and continue
        fig_name = ['corr_TEP-SICIxMEP-SICI_ranked_' varnames{row(a)} '_' varnames{col(a)}];
        savefig([folder_figures '\' fig_name '.fig'])
        saveas(fig, [folder_figures '\' fig_name  '.png'])  
        figure_counter = figure_counter + 1;
    end
end
clear peaks_SICI varnames m p data_cor alpha_cor cor_coef cor_p row col data_model fig plot_cor fig_name T text_pos a temp 

%% ***
% M1: hostogram - peak amplitude
fig = figure(figure_counter);
for p = 1:length(peaks_M1)
    for m = 1:length(medication)        
        % launch the plot
        subplot(length(medication), length(peaks_M1), (m-1)*length(peaks_M1) + p)
        hold on
        
        % plot the histogram with a normal distribution fit
        h = histfit(data_M1(:, m, p, 1), 10, 'kernel');
        h(1).FaceColor = colours((m-1)*2 + 1, :); h(1).EdgeColor = [1 1 1];
        h(2).Color = [0 0 0]; h(2).LineWidth = 1;
        
        % get limits
        yl = get(gca, 'ylim'); xl = get(gca, 'xlim');
        
        % add mean of the distribution
        mu = mean(data_M1(:, m, p, 1));
        text(xl(1), yl(2)*-0.1, sprintf('mu = %1.3f', mu), 'FontSize', 10, 'FontWeight', 'bold')
        
        % add skewness and kurtosis
        s = skewness(data_M1(:, m, p, 1)); k = kurtosis(data_M1(:, m, p, 1));
        text(xl(1), yl(2)*-0.2, sprintf('%1.3f, %1.3f', s, k), 'FontSize', 10)
            
        % set annotation
        if p == 1
            text(xl(1) - (xl(2) - xl(1)), yl(2)*0.7, medication{m}, 'FontSize', 10)
            clear yl xl
        end
                       
        % set title
        if m == 1
            title(peaks_M1{p})
        end
    end
end   
clear p m h yl xl mu s k

% parameters
fig.Position = [500 400 750 500];
sgtitle('M1 - change in peak amplitude')

% update counter
figure_counter = figure_counter + 1;

% M1: QQ plot - peak amplitude
fig = figure(figure_counter);
for p = 1:length(peaks_M1)
    for m = 1:length(medication)        
        % launch the plot
        subplot(length(medication), length(peaks_M1), (m-1)*length(peaks_M1) + p)
        hold on
        
        % plot the histogram with a normal distribution fit
        q = qqplot(data_M1(:, m, p, 1));
        title(''); xlabel(''); ylabel('')
        
        % get limits
        yl = get(gca, 'ylim'); xl = get(gca, 'xlim');
      
        % set annotation
        if p == 1
            text(xl(1) - (xl(2) - xl(1))*0.7, yl(2)*0.7, medication{m}, 'FontSize', 10)
            clear yl xl
        end
                       
        % set title
        if m == 1
            title(peaks_M1{p})
        end
    end
end   
clear p m q yl xl

% parameters
fig.Position = [200 400 1000 500];
sgtitle('M1 - change in peak amplitude')

% update counter
figure_counter = figure_counter + 1;

% AG: histogram - peak amplitude
fig = figure(figure_counter);
for p = 1:length(peaks_AG)
    for m = 1:length(medication)        
        % launch the plot
        subplot(length(medication), length(peaks_AG), (m-1)*length(peaks_AG) + p)
        hold on
        
        % plot the histogram with a normal distribution fit
        h = histfit(data_AG(:, m, p, 1), 10, 'kernel');
        h(1).FaceColor = colours((m-1)*2 + 1, :); h(1).EdgeColor = [1 1 1];
        h(2).Color = [0 0 0]; h(2).LineWidth = 1;
        
        % get limits
        yl = get(gca, 'ylim'); xl = get(gca, 'xlim');
        
        % add mean of the distribution
        mu = mean(data_AG(:, m, p, 1));
        text(xl(1), yl(2)*-0.1, sprintf('mu = %1.3f', mu), 'FontSize', 10, 'FontWeight', 'bold')
        
        % add skewness and kurtosis
        s = skewness(data_AG(:, m, p, 1)); k = kurtosis(data_AG(:, m, p, 1));
        text(xl(1), yl(2)*-0.2, sprintf('%1.3f, %1.3f', s, k), 'FontSize', 10)
            
        % set annotation
        if p == 1
            text(xl(1) - (xl(2) - xl(1)), yl(2)*0.7, medication{m}, 'FontSize', 10)
            clear yl xl
        end
                       
        % set title
        if m == 1
            title(peaks_AG{p})
        end
    end
end   
clear p m h yl xl mu s k

% parameters
fig.Position = [500 400 750 500];
sgtitle('AG - change in peak amplitude')

% update counter
figure_counter = figure_counter + 1;

% AG: QQ plot - peak amplitude
fig = figure(figure_counter);
for p = 1:length(peaks_AG)
    for m = 1:length(medication)        
        % launch the plot
        subplot(length(medication), length(peaks_AG), (m-1)*length(peaks_AG) + p)
        hold on
        
        % plot the histogram with a normal distribution fit
        q = qqplot(data_AG(:, m, p, 1));
        title(''); xlabel(''); ylabel('')
        
        % get limits
        yl = get(gca, 'ylim'); xl = get(gca, 'xlim');
      
        % set annotation
        if p == 1
            text(xl(1) - (xl(2) - xl(1))*0.7, yl(2)*0.7, medication{m}, 'FontSize', 10)
            clear yl xl
        end
                       
        % set title
        if m == 1
            title(peaks_AG{p})
        end
    end
end   
clear p m q yl xl 

% parameters
fig.Position = [200 400 1000 500];
sgtitle('AG - change in peak amplitude')

% update counter
figure_counter = figure_counter + 1;
clear fig




