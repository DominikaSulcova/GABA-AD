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

% average TEPs
load([folder_results '\GABA_' group '_variables\GABA_' group '_M1_TEPs.mat'], 'GABA_data_mean')
TEP_M1 = GABA_data_mean;
disp(['M1 TEPs - datasize: ' num2str(size(TEP_M1))])

load([folder_results '\GABA_' group '_variables\GABA_' group '_AG_TEPs.mat'], 'GABA_data_mean')
TEP_AG = GABA_data_mean;
disp(['AG TEPs - datasize: ' num2str(size(TEP_AG))])

TEP_mean(1, :, :, :, :) = squeeze(TEP_M1(:, :, 1, 1:30, :));
TEP_mean(2, :, :, :, :) = squeeze(TEP_M1(:, :, 2, 1:30, :));
TEP_mean(3, :, :, :, :) = squeeze(TEP_AG(:, :, 1:30, :));
TEP_mean(4, :, :, :, :) = squeeze(TEP_M1(:, :, 3, 1:30, :));
clear GABA_data_mean TEP_AG TEP_M1

% SICI
load([folder_results '\GABA_' group '_variables\GABA_' group '_M1_TEPs.mat'], 'GABA_SICI')

%% ) BASELINE TEPs
% time axis
x = [-50:0.5:300];

% plot the baseline butterfly plots for all stimuli
for s = 1:3
    % prepare baseline data data (placebo session)
    data_visual = squeeze(TEP_mean(s, 1, 1, :, :));

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

%% ) BASELINE SICI - AMPLITUDE CHANGE
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

%% ) BASELINE SICI - TEP-MEP CORRELATION
% ----- decide output parameters -----
peaks_SICI = {'N17' 'P60' 'N100'}; 
% ------------------------------------
% get the data
for m = 1:length(medication)
    for p = 1:length(participant)    
        % TEP SICI = x 
        x_corr(:, (m-1)*length(participant) + p) = GABA_YC_results.SICI(m).TEP.pre(p, contains(peaks_M1, peaks_SICI))';
        
        % MEP SICI = y
        y_corr((m-1)*length(participant) + p) = 100 - GABA_YC_results.SICI(m).MEP.pre(p);
        
        % choose colours
        marker_col((m-1)*length(participant) + p, :) = colours((m-1)*2 + 2, :);
    end
end
clear m p

% calculate and plot for each component of interest
for p = 1:length(peaks_SICI)
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

    % plot data + regression line
    fig = figure(figure_counter);
    hold on
    plot_corr(data_model_ranked, data_corr_ranked, marker_col, 'Spearman')

    % save the figure   
    figure_name = ['corr_SICI_MEPx' peaks_SICI{p}];
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.png'])
    
    % update figure counter
    figure_counter = figure_counter + 1;
end
clear x_corr y_corr marker_col data_corr data_corr_ranked data_model_ranked fig figure_name

for p = 1:length(peaks_SICI)
        % fit
    [curvefit,gof,output] = fit(x, y, 'exp1');
    
    % plot
    figure(figure_counter)
    subplot(2, 2, 1)
    plot(curvefit,x,y)
    
    subplot(2, 2, 2)
    scatter(x,log(y))  
    
    subplot(2, 2, 3)
    plot(curvefit,x,y,'Residuals')
    
    subplot(2, 2, 4)
    plot(curvefit,x,y,'predfunc')
    
    % compute coefficients
    n = length(x);
    y2 = log(y);
    j = sum(x);
    k = sum(y2);
    r2 = sum(x .* y2);
    m = sum(y2.^2);
    c = f.b * (r2 - j * k / n);
    d = m - k^2 / n;
    corr = sqrt(c/d);
    std_err = sqrt((d - c) / (n - 2)); 
    
end
clear x_corr y_corr marker_col p x y curvefit n y2 r2 j k m c d corr std_err gof output

%% ) BASELINE SICI - SUBTRACTION
% ----- decide output parameters -----
electrode = {'C3'}; 
% ------------------------------------
% load default labels
load([folder_git '\GABA_header_default.mat'])
labels = {header.chanlocs.labels};
clear header

% time axis
x = [-50:0.5:300];

% plot SICI for each chosen electrode
for e = 1:length(electrode)
    % prepare baseline data data at C3 (placebo session)
    e_n = find(strcmp(labels, electrode{e}));
    data_visual = squeeze(GABA_SICI.mean(1, 1, e_n, :))';
    CI_visual = squeeze(GABA_SICI.CI(1, 1, e_n, :))';

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
    ylabel('amplitude (\muV)')
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
    saveas(fig, [folder_figures '\' figure_name '.png'])
end

% update figure counter
figure_counter = figure_counter + 1 ;

clear electrode e e_n x data_visual CI_visual fig yl P F lgd figure_name 

%% FUNCTIONS
function plot_corr(data_model, data_corr, marker_col, corr_type)
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

% add annotations
text_pos = [0.90 0.78];
rectangle('Position', [2, 27, 18, 12], 'FaceColor',[1 1 1], 'EdgeColor', [0.5 0.5 0.5])
T(1) = text(0.1, text_pos(1), sprintf( 'y = %1.3f * x', data_model.Coefficients.Estimate(2)), 'Units', 'Normalized');
T(2) = text(0.1, text_pos(2), sprintf('R^2 = %1.3f', data_model.Rsquared.Adjusted), 'Units', 'Normalized');
set(T(1), 'fontsize', 20, 'fontangle', 'italic', 'fontweight', 'bold'); 
set(T(2), 'fontsize', 20); 
end


