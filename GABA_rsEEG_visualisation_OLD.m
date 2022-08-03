%% RS-EEG - FREQUENCY CONTENT VISUALISATION - relative amplitude 
% Written by Dominika for GABA-AD project (2021)
% 
% Plots mean amplitude values separately for higer frequencies [4 - 45]Hz
% and lower = delta frequencies [0.1 - 4]Hz  
% 1) Extracts mean values
% 2) Plots baseline freqency content for all ROIs
% 3) Plots baseline against post-medication freqency content for all ROIs
%       + some close-ups


%% parameters
clear all; clc;

% dataset
participant = 1:20;  
medication = {'placebo' 'alprazolam'};      
time = {'pre' 'post'};
condition = {'open' 'closed'};
prefix_1 = 'avg fft han icfilt assigned ica';
prefix_2 = 'crop_120 45Hz_but4 but reref chan-select ds YC';

% ROIs
ROI = struct;
ROI(1).area = 'frontal'; ROI(2).area = 'central'; ROI(3).area = 'left_temporal'; ROI(4).area = 'right_temporal'; ROI(5).area = 'occipital'; 
ROI(1).electrodes = {'Fp1' 'Fp2' 'Fz' 'F3' 'F4' 'F7' 'F8'};
ROI(2).electrodes = {'FC1' 'FC2' 'Cz' 'C1' 'C2' 'CP1' 'CP2'};
ROI(3).electrodes = {'FC5' 'T7' 'C3' 'CP5'};
ROI(4).electrodes = {'FC6' 'T8' 'C4' 'CP6'};
ROI(5).electrodes = {'P3' 'P4' 'Pz' 'P7' 'P8' 'O1' 'O2' 'Iz'};
load('labels.mat')
labels = labels_TS(1:30);
clear labels_CS labels_TS

% TOIs
TOI = struct;
TOI(1).band = 'delta'; TOI(2).band = 'theta'; TOI(3).band = 'alpha'; TOI(4).band = 'beta'; TOI(5).band = 'gamma';
TOI(1).window = [0.1, 4]; TOI(2).window = [4, 8]; TOI(3).window = [8, 13]; TOI(4).window = [13, 30]; TOI(5).window = [30, 45];
TOI(1).sign = '\delta'; TOI(2).sign = '\theta'; TOI(3).sign = '\alpha'; TOI(4).sign = '\beta'; TOI(5).sign = '\gamma';

% import header basis
load([prefix_1 ' high ' prefix_2 '1 placebo pre EEGcont open.lw6'], '-mat')
header_high = header;
load([prefix_1 ' low ' prefix_2 '1 placebo pre EEGcont open.lw6'], '-mat')
header_low = header;
clear header

% graphic parameters
figure_counter = 1;
alpha = 0.2;

% statistics
z = 1.96;

%% load the data
load('rsEEG_data_high.mat')
load('rsEEG_data_low.mat') 

%% HIGHER FREQUENCIES: 1) calculate mean values
% preapre outcome variables
mean_high = [];
CI_high = [];
mean_bl_high = [];
CI_bl_high = [];

% mean values for each condition separately
for m = 1:size(data_high, 1)
    for t = 1:size(data_high, 2)
        for c = 1:size(data_high, 3)
            for r = 1:size(data_high, 5)
                for i = 1:size(data_high, 6)
                    mean_high(m, t, c, r, i) = mean(data_high(m, t, c, :, r, i));
                    CI_high(m, t, c, r, i) = (std(data_high(m, t, c, :, r, i))/sqrt(size(data_high, 4)));%*z;
                end                
            end
        end
    end
end
datasize_mean = [size(mean_high, 1) size(mean_high, 2) size(mean_high, 3) size(mean_high, 4) size(mean_high, 5)];

% mean values for baseline conditions together 
for c = 1:datasize(3)
    for r = 1:datasize(5)
        for i = 1:datasize(6)
            data_bl = [squeeze(data_high(1, 1, c, :, r, i)); squeeze(data_high(2, 1, c, :, r, i))];
            mean_bl_high(c, r, i) = mean(data_bl);
            CI_bl_high(c, r, i) = (std(data_bl)/sqrt(length(data_bl)))*z;
        end                
    end
end

%% HIGHER FREQUENCIES: 2) plot baseline
% decide colour scheme
answer = questdlg('Would you like to set a new colour scheme?', 'Visualization - colours',...
    'Yes', 'No, load the preset', 'No, load the preset');
switch answer
    case 'Yes'
        for CS = 1:numel(ROI)
            colours(CS, :) = uisetcolor(['Choose colour for the ' ROI(CS).area ' ROI:'])
        end
    case 'No, load the preset'
        load('colours.mat');
end
clear answer

% determine x axis
x_start = header_high.xstart;
x_end = header_high.datasize(6) * header_high.xstep + header_high.xstart - header_high.xstep;
x_high = [x_start : header_high.xstep : x_end];

% ----- plot signal with open X closed eyes for each ROI -----
for r = 1:datasize(5)
    % choose data
    data_visual_open = squeeze(mean_bl_high(1, r, :))'; data_visual_closed = squeeze(mean_bl_high(2, r, :))';  
    CI_visual_open = squeeze(CI_bl_high(1, r, :))'; CI_visual_closed = squeeze(CI_bl_high(2, r, :))'; 
    
    % set limits of the figure
    fig = figure(figure_counter);
    plot(x_high, data_visual_closed + CI_visual_closed, 'b:', 'LineWidth', 0.5)
    yl = get(gca, 'ylim');
    clf
    ylim([-0.2, yl(2)])
    xlim([x_start, x_end])
    hold on

    % plot lines between frequency bands (= TOIs)
    h1 = line([TOI(2).window(2), TOI(2).window(2)], [-0.2, yl(2)], 'Color', [0.9 0.9 0.9], 'LineWidth', 0.8); 
    h2 = line([TOI(3).window(2), TOI(3).window(2)], [-0.2, yl(2)], 'Color', [0.9 0.9 0.9], 'LineWidth', 0.8);
    h3 = line([TOI(4).window(2), TOI(4).window(2)], [-0.2, yl(2)], 'Color', [0.9 0.9 0.9], 'LineWidth', 0.8);
    h4 = line([x_start, x_end], [0, 0], 'Color', [0, 0, 0], 'LineWidth', 0.5, 'LineStyle', ':');
    hold on
    
    % add names of frequency bands
    for f = 1:4
        text((TOI(f + 1).window(1) + TOI(f + 1).window(2))/2, -0.09, TOI(f + 1).sign, 'fontsize', 16, 'color', [0.6 0.6 0.6])
        hold on
    end
    
    % plot data - eyes open   
    p1 = plot(x_high,  data_visual_open, 'Color', colours(r, :), 'LineWidth', 2.5, 'LineStyle', ':');
    
    xlabel('frequency (Hz)')
    ylabel('relative amplitude (µV)')
    set(gca, 'FontSize', 16)
    title([ROI(r).area ' region'], 'FontWeight', 'bold', 'FontSize', 18)
    hold on

    % add CI
    f1 = fill([x_high fliplr(x_high)],[data_visual_open + CI_visual_open fliplr(data_visual_open - CI_visual_open)], ...
        colours(r, :), 'FaceAlpha', alpha, 'linestyle', 'none');
    set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
    hold on

    % plot data - eyes closed   
    p2 = plot(x_high,  data_visual_closed, 'Color', colours(r, :), 'LineWidth', 2.5);
    hold on

    % add CI
    f2 = fill([x_high fliplr(x_high)],[data_visual_closed + CI_visual_closed fliplr(data_visual_closed - CI_visual_closed)], ...
        colours(r, :), 'FaceAlpha', alpha, 'linestyle', 'none');
    set(get(get(f2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
    
    % add legend 
    legend([p1, p2], {'eyes open' 'eyes closed'}, 'FontSize', 16, 'position', [0.55, 0.55, 0.3 0.2], 'edgecolor', [0.6, 0.6, 0.6])
    hold off

    % save figure
    figure_name = ['rsEEG_baseline_' ROI(r).area];
    savefig(figure_name)
    saveas(fig, [figure_name '.png'])
    pause(3)

    % update counter
    figure_counter = figure_counter + 1;
end

% ----- plot signal with closed eyes for all ROIs together -----
% choose data
data_visual_closed = squeeze(mean_bl_high(2, :, :));  

% launch the figure
fig = figure(figure_counter);
xlim([x_start, x_end])
ylim([-0.2, 1.2])
% r1 = rectangle('Position', [TOI(3).window(1), -0.2, TOI(3).window(2) - TOI(3).window(1), 1.4], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
h1 = line([TOI(2).window(2), TOI(2).window(2)], [-0.2, 1.2], 'Color', [0.9 0.9 0.9], 'LineWidth', 0.8); 
h1 = line([TOI(3).window(2), TOI(3).window(2)], [-0.2, 1.2], 'Color', [0.9 0.9 0.9], 'LineWidth', 0.8); 
h1 = line([TOI(4).window(2), TOI(4).window(2)], [-0.2, 1.2], 'Color', [0.9 0.9 0.9], 'LineWidth', 0.8); 
h4 = line([x_start, x_end], [0, 0], 'Color', [0, 0, 0], 'LineWidth', 0.5, 'LineStyle', ':');
hold on

% plot the data
for r = 1:size(data_visual_closed, 1)
    pl(r) = plot(x_high, data_visual_closed(r, :), 'Color', colours(r, :), 'LineWidth', 2);
    hold on
end

% add names of frequency bands
for f = 1:4
    text((TOI(f + 1).window(1) + TOI(f + 1).window(2))/2, -0.09, TOI(f + 1).sign, 'fontsize', 16, 'color', [0.6 0.6 0.6])
    hold on
end

% add parameters  
xlabel('frequency (Hz)')
ylabel('relative amplitude (µV)')
set(gca, 'FontSize', 16)
title('eyes closed', 'FontWeight', 'bold', 'FontSize', 18)
hold on

% add legend 
legend(pl(1:5), {ROI(1:5).area}, 'FontSize', 16, 'position', [0.55, 0.55, 0.3, 0.3], 'edgecolor', [0.6 0.6 0.6])
hold off

% save figure
figure_name = ['rsEEG_baseline_closed'];
savefig(figure_name)
saveas(fig, [figure_name '.png'])
pause(3)

% update counter
figure_counter = figure_counter + 1;

% ----- plot signal with open eyes for all ROIs together -----
% choose data
data_visual_open = squeeze(mean_bl_high(1, :, :));  

% launch the figure
fig = figure(figure_counter);
xlim([x_start, x_end])
ylim([-0.2, 1.2])
% r1 = rectangle('Position', [TOI(3).window(1), -0.2, TOI(3).window(2) - TOI(3).window(1), 1.4], 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
h1 = line([TOI(2).window(2), TOI(2).window(2)], [-0.2, 1.2], 'Color', [0.9 0.9 0.9], 'LineWidth', 0.8); 
h1 = line([TOI(3).window(2), TOI(3).window(2)], [-0.2, 1.2], 'Color', [0.9 0.9 0.9], 'LineWidth', 0.8); 
h1 = line([TOI(4).window(2), TOI(4).window(2)], [-0.2, 1.2], 'Color', [0.9 0.9 0.9], 'LineWidth', 0.8); 
h4 = line([x_start, x_end], [0, 0], 'Color', [0, 0, 0], 'LineWidth', 0.5, 'LineStyle', ':');
hold on

% plot the data
for r = 1:size(data_visual_open, 1)
    pl(r) = plot(x_high, data_visual_open(r, :), 'Color', colours(r, :), 'LineWidth', 2);
    hold on
end

% add names of frequency bands
for f = 1:4
    text((TOI(f + 1).window(1) + TOI(f + 1).window(2))/2, -0.09, TOI(f + 1).sign, 'fontsize', 16, 'color', [0.6 0.6 0.6])
    hold on
end

% add parameters  
xlabel('frequency (Hz)')
ylabel('relative amplitude (µV)')
set(gca, 'FontSize', 16)
title('eyes open', 'FontWeight', 'bold', 'FontSize', 18)
hold on

% add legend 
legend(pl(1:5), {ROI(1:5).area}, 'FontSize', 16, 'position', [0.55, 0.55, 0.3, 0.3], 'edgecolor', [0.6 0.6 0.6])
hold off

% save figure
figure_name = ['rsEEG_baseline_open'];
savefig(figure_name)
saveas(fig, [figure_name '.png'])
pause(3)

% update counter
figure_counter = figure_counter + 1;

%% HIGHER FREQUENCIES: 3) plot difference pre X post
% decide colour scheme
answer = questdlg('Would you like to set a new colour scheme?', 'Visualization - colours',...
    'Yes', 'No, load the preset', 'No, load the preset');
switch answer
    case 'Yes'
        for CS = 1:4
            colours2(CS, :) = uisetcolor(['Choose colour ' num2str(CS)])
        end
    case 'No, load the preset'
        load('colours2.mat');
end
clear answer


% ----- plot baseline signal for both medications -----
% choose data
data_pla_visual_open = squeeze(mean_high(2, 1, 1, :, :)); data_pla_visual_closed = squeeze(mean_high(1, 1, 2, :, :));  
CI_pla_visual_open = squeeze(CI_high(2, 1, 1, :, :)); CI_pla_visual_closed = squeeze(CI_high(1, 1, 2, :, :)); 
data_alp_visual_open = squeeze(mean_high(1, 1, 1, :, :)); data_alp_visual_closed = squeeze(mean_high(2, 1, 2, :, :));  
CI_alp_visual_open = squeeze(CI_high(1, 1, 1, :, :)); CI_alp_visual_closed = squeeze(CI_high(2, 1, 2, :, :)); 

% plot occipital region
for r = 5  
    % set limits of the figure
    fig = figure(figure_counter);
    plot(x_high, data_pla_visual_closed(r, :), 'b:', 'LineWidth', 0.5)
    yl = get(gca, 'ylim');
    clf
    ylim([-0.2, yl(2)])
    xlim([x_start, x_end])
    hold on

    % plot lines between frequency bands (= TOIs)
    h1 = line([TOI(2).window(2), TOI(2).window(2)], [-0.2, yl(2)], 'Color', [0.9 0.9 0.9], 'LineWidth', 0.8); 
    h2 = line([TOI(3).window(2), TOI(3).window(2)], [-0.2, yl(2)], 'Color', [0.9 0.9 0.9], 'LineWidth', 0.8);
    h3 = line([TOI(4).window(2), TOI(4).window(2)], [-0.2, yl(2)], 'Color', [0.9 0.9 0.9], 'LineWidth', 0.8);
    h4 = line([x_start, x_end], [0, 0], 'Color', [0, 0, 0], 'LineWidth', 0.5, 'LineStyle', ':');
    hold on
    
    % add names of frequency bands
    for f = 1:4
        text((TOI(f + 1).window(1) + TOI(f + 1).window(2))/2, -0.09, TOI(f + 1).sign, 'fontsize', 16, 'color', [0.6 0.6 0.6])
        hold on
    end
    
    % plot data & CI - placebo, eyes open     
    p1 = plot(x_high,  data_pla_visual_open(r, :), 'Color', colours2(1, :), 'LineWidth', 2.5, 'LineStyle', ':');    
    xlabel('frequency (Hz)')
    ylabel('relative amplitude (µV)')
    set(gca, 'FontSize', 16)
    title([ROI(r).area ' region'], 'FontWeight', 'bold', 'FontSize', 18)
    hold on
%     f1 = fill([x_high fliplr(x_high)],[data_pla_visual_open(r, :) + CI_pla_visual_open(r, :) fliplr(data_pla_visual_open(r, :) - CI_pla_visual_open(r, :))], ...
%         colours2(1, :), 'FaceAlpha', alpha, 'linestyle', 'none');
%     set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
%     hold on

    % plot data & CI - placebo, eyes closed     
    p2 = plot(x_high,  data_pla_visual_closed(r, :), 'Color', colours2(1, :), 'LineWidth', 2.5);    
    hold on
%     f2 = fill([x_high fliplr(x_high)],[data_pla_visual_closed(r, :) + CI_pla_visual_closed(r, :) fliplr(data_pla_visual_closed(r, :) - CI_pla_visual_closed(r, :))], ...
%         colours2(1, :), 'FaceAlpha', alpha, 'linestyle', 'none');
%     set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
%     hold on
    
    % plot data & CI - alprazolam, eyes open   
    p3 = plot(x_high,  data_alp_visual_open(r, :), 'Color', colours2(3, :), 'LineWidth', 2.5, 'linestyle', ':');    
    hold on
%     f3 = fill([x_high fliplr(x_high)],[data_alp_visual_open(r, :) + CI_alp_visual_open(r, :) fliplr(data_alp_visual_open(r, :) - CI_alp_visual_open(r, :))], ...
%         colours2(3, :), 'FaceAlpha', alpha, 'linestyle', 'none');
%     set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
%     hold on
    
    % plot data & CI - alprazolam, eyes closed   
    p4 = plot(x_high,  data_alp_visual_closed(r, :), 'Color', colours2(3, :), 'LineWidth', 2.5);    
    hold on
%     f4 = fill([x_high fliplr(x_high)],[data_alp_visual_closed(r, :) + CI_alp_visual_closed(r, :) fliplr(data_alp_visual_closed(r, :) - CI_alp_visual_closed(r, :))], ...
%         colours2(3, :), 'FaceAlpha', alpha, 'linestyle', 'none');
%     set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
%     hold on
    
    % add legend 
    legend([p1, p2, p3, p4], {'placebo - eyes open' 'placebo - eyes closed' 'alprazolam - eyes open' 'alprazolam - eyes closed'}, ...
        'FontSize', 16, 'position', [0.48, 0.55, 0.3 0.2], 'edgecolor', [0.6, 0.6, 0.6])
    hold off

    % save figure
    figure_name = ['rsEEG_baseline_both_' ROI(r).area];
    savefig(figure_name)
    saveas(fig, [figure_name '.png'])
    pause(3)

    % update counter
    figure_counter = figure_counter + 1;
end


% ----- plot pre X post medication for placebo -----
% choose data
data_pre_visual_open = squeeze(mean_high(1, 1, 1, :, :)); data_pre_visual_closed = squeeze(mean_high(1, 1, 2, :, :));  
CI_pre_visual_open = squeeze(CI_high(1, 1, 1, :, :)); CI_pre_visual_closed = squeeze(CI_high(1, 1, 2, :, :)); 
data_post_visual_open = squeeze(mean_high(1, 2, 1, :, :)); data_post_visual_closed = squeeze(mean_high(1, 2, 2, :, :));  
CI_post_visual_open = squeeze(CI_high(1, 2, 1, :, :)); CI_post_visual_closed = squeeze(CI_high(1, 2, 2, :, :)); 

% plot all regions
for r = 1:datasize(5)
    % set limits of the figure
    fig = figure(figure_counter);
    plot(x_high, data_pre_visual_closed(r, :) + CI_pre_visual_closed(r, :), 'b:', 'LineWidth', 0.5)
    yl = get(gca, 'ylim');
    clf
    ylim([-0.2, yl(2)])
    xlim([x_start, x_end])
    hold on

    % plot lines between frequency bands (= TOIs)
    h1 = line([TOI(2).window(2), TOI(2).window(2)], [-0.2, yl(2)], 'Color', [0.9 0.9 0.9], 'LineWidth', 0.8); 
    h2 = line([TOI(3).window(2), TOI(3).window(2)], [-0.2, yl(2)], 'Color', [0.9 0.9 0.9], 'LineWidth', 0.8);
    h3 = line([TOI(4).window(2), TOI(4).window(2)], [-0.2, yl(2)], 'Color', [0.9 0.9 0.9], 'LineWidth', 0.8);
    h4 = line([x_start, x_end], [0, 0], 'Color', [0, 0, 0], 'LineWidth', 0.5, 'LineStyle', ':');
    hold on
    
    % add names of frequency bands
    for f = 1:4
        text((TOI(f + 1).window(1) + TOI(f + 1).window(2))/2, -0.09, TOI(f + 1).sign, 'fontsize', 16, 'color', [0.6 0.6 0.6])
        hold on
    end
    
    % plot data & CI - pre, eyes open     
    p1 = plot(x_high,  data_pre_visual_open(r, :), 'Color', colours2(1, :), 'LineWidth', 2.5, 'LineStyle', ':');    
    xlabel('frequency (Hz)')
    ylabel('relative amplitude (µV)')
    set(gca, 'FontSize', 16)
    title(['placebo - ' ROI(r).area ' region'], 'FontWeight', 'bold', 'FontSize', 18)
    hold on
    f1 = fill([x_high fliplr(x_high)],[data_pre_visual_open(r, :) + CI_pre_visual_open(r, :) fliplr(data_pre_visual_open(r, :) - CI_pre_visual_open(r, :))], ...
        colours2(1, :), 'FaceAlpha', alpha, 'linestyle', 'none');
    set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
    hold on

    % plot data & CI - pre, eyes closed     
    p2 = plot(x_high,  data_pre_visual_closed(r, :), 'Color', colours2(1, :), 'LineWidth', 2.5);    
    hold on
    f2 = fill([x_high fliplr(x_high)],[data_pre_visual_closed(r, :) + CI_pre_visual_closed(r, :) fliplr(data_pre_visual_closed(r, :) - CI_pre_visual_closed(r, :))], ...
        colours2(1, :), 'FaceAlpha', alpha, 'linestyle', 'none');
    set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
    hold on
    
    % plot data & CI - post, eyes open   
    p3 = plot(x_high,  data_post_visual_open(r, :), 'Color', colours2(2, :), 'LineWidth', 2.5, 'linestyle', ':');    
    hold on
    f3 = fill([x_high fliplr(x_high)],[data_post_visual_open(r, :) + CI_post_visual_open(r, :) fliplr(data_post_visual_open(r, :) - CI_post_visual_open(r, :))], ...
        colours2(2, :), 'FaceAlpha', alpha, 'linestyle', 'none');
    set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
    hold on
    
    % plot data & CI - alprazolam, eyes closed   
    p4 = plot(x_high,  data_post_visual_closed(r, :), 'Color', colours2(2, :), 'LineWidth', 2.5);    
    hold on
    f4 = fill([x_high fliplr(x_high)],[data_post_visual_closed(r, :) + CI_post_visual_closed(r, :) fliplr(data_post_visual_closed(r, :) - CI_post_visual_closed(r, :))], ...
        colours2(2, :), 'FaceAlpha', alpha, 'linestyle', 'none');
    set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
    hold on
    
    % add legend 
    legend([p1, p2, p3, p4], {'baseline - eyes open' 'baseline - eyes closed' 'post - eyes open' 'post - eyes closed'}, ...
        'FontSize', 16, 'position', [0.48, 0.55, 0.3 0.3], 'edgecolor', [0.6, 0.6, 0.6])
    hold off

    % save figure
    figure_name = ['rsEEG_placebo_' ROI(r).area];
    savefig(figure_name)
    saveas(fig, [figure_name '.png'])
    pause(3)

    % update counter
    figure_counter = figure_counter + 1;
end

% ----- plot pre X post medication for alprazolam -----
% choose data
data_pre_visual_open = squeeze(mean_high(2, 1, 1, :, :)); data_pre_visual_closed = squeeze(mean_high(2, 1, 2, :, :));  
CI_pre_visual_open = squeeze(CI_high(2, 1, 1, :, :)); CI_pre_visual_closed = squeeze(CI_high(2, 1, 2, :, :)); 
data_post_visual_open = squeeze(mean_high(2, 2, 1, :, :)); data_post_visual_closed = squeeze(mean_high(2, 2, 2, :, :));  
CI_post_visual_open = squeeze(CI_high(2, 2, 1, :, :)); CI_post_visual_closed = squeeze(CI_high(2, 2, 2, :, :)); 

% plot all regions
for r = 1:datasize(5)
    % set limits of the figure
    fig = figure(figure_counter);
    plot(x_high, data_pre_visual_closed(r, :) + CI_pre_visual_closed(r, :), 'b:', 'LineWidth', 0.5)
    yl = get(gca, 'ylim');
    clf
    ylim([-0.2, yl(2)])
    xlim([x_start, x_end])
    hold on

    % plot lines between frequency bands (= TOIs)
    h1 = line([TOI(2).window(2), TOI(2).window(2)], [-0.2, yl(2)], 'Color', [0.9 0.9 0.9], 'LineWidth', 0.8); 
    h2 = line([TOI(3).window(2), TOI(3).window(2)], [-0.2, yl(2)], 'Color', [0.9 0.9 0.9], 'LineWidth', 0.8);
    h3 = line([TOI(4).window(2), TOI(4).window(2)], [-0.2, yl(2)], 'Color', [0.9 0.9 0.9], 'LineWidth', 0.8);
    h4 = line([x_start, x_end], [0, 0], 'Color', [0, 0, 0], 'LineWidth', 0.5, 'LineStyle', ':');
    hold on
    
    % add names of frequency bands
    for f = 1:4
        text((TOI(f + 1).window(1) + TOI(f + 1).window(2))/2, -0.09, TOI(f + 1).sign, 'fontsize', 16, 'color', [0.6 0.6 0.6])
        hold on
    end
    
    % plot data & CI - pre, eyes open     
    p1 = plot(x_high,  data_pre_visual_open(r, :), 'Color', colours2(3, :), 'LineWidth', 2.5, 'LineStyle', ':');    
    xlabel('frequency (Hz)')
    ylabel('relative amplitude (µV)')
    set(gca, 'FontSize', 16)
    title(['alprazolam - ' ROI(r).area ' region'], 'FontWeight', 'bold', 'FontSize', 18)
    hold on
    f1 = fill([x_high fliplr(x_high)],[data_pre_visual_open(r, :) + CI_pre_visual_open(r, :) fliplr(data_pre_visual_open(r, :) - CI_pre_visual_open(r, :))], ...
        colours2(3, :), 'FaceAlpha', alpha, 'linestyle', 'none');
    set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
    hold on

    % plot data & CI - pre, eyes closed     
    p2 = plot(x_high,  data_pre_visual_closed(r, :), 'Color', colours2(3, :), 'LineWidth', 2.5);    
    hold on
    f2 = fill([x_high fliplr(x_high)],[data_pre_visual_closed(r, :) + CI_pre_visual_closed(r, :) fliplr(data_pre_visual_closed(r, :) - CI_pre_visual_closed(r, :))], ...
        colours2(3, :), 'FaceAlpha', alpha, 'linestyle', 'none');
    set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
    hold on
    
    % plot data & CI - post, eyes open   
    p3 = plot(x_high,  data_post_visual_open(r, :), 'Color', colours2(4, :), 'LineWidth', 2.5, 'linestyle', ':');    
    hold on
    f3 = fill([x_high fliplr(x_high)],[data_post_visual_open(r, :) + CI_post_visual_open(r, :) fliplr(data_post_visual_open(r, :) - CI_post_visual_open(r, :))], ...
        colours2(4, :), 'FaceAlpha', alpha, 'linestyle', 'none');
    set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
    hold on
    
    % plot data & CI - alprazolam, eyes closed   
    p4 = plot(x_high,  data_post_visual_closed(r, :), 'Color', colours2(4, :), 'LineWidth', 2.5);    
    hold on
    f4 = fill([x_high fliplr(x_high)],[data_post_visual_closed(r, :) + CI_post_visual_closed(r, :) fliplr(data_post_visual_closed(r, :) - CI_post_visual_closed(r, :))], ...
        colours2(4, :), 'FaceAlpha', alpha, 'linestyle', 'none');
    set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
    hold on
    
    % add legend 
    legend([p1, p2, p3, p4], {'baseline - eyes open' 'baseline - eyes closed' 'post - eyes open' 'post - eyes closed'}, ...
        'FontSize', 16, 'position', [0.48, 0.55, 0.3 0.3], 'edgecolor', [0.6, 0.6, 0.6])
    hold off

    % save figure
    figure_name = ['rsEEG_alprazolam_' ROI(r).area];
    savefig(figure_name)
    saveas(fig, [figure_name '.png'])
    pause(3)

    % update counter
    figure_counter = figure_counter + 1;
end

% plot a close up - increase of eyes-open alpha over the frontal region
% choose data
data_closeup_pre = data_pre_visual_open(1, 2/header_high.xstep:11/header_high.xstep); data_closeup_post = data_post_visual_open(1, 2/header_high.xstep:11/header_high.xstep);
CI_closeup_pre = CI_pre_visual_open(1, 2/header_high.xstep:11/header_high.xstep); CI_closeup_post = CI_post_visual_open(1, 2/header_high.xstep:11/header_high.xstep);

% adapt x
x_closeup = x_high(2/header_high.xstep:11/header_high.xstep);

% set limits of the figure
fig = figure(figure_counter);
plot(x_closeup, data_closeup_post + CI_closeup_post, 'b:', 'LineWidth', 0.5)
yl = get(gca, 'ylim');
clf
ylim([0, yl(2)])
xlim([6, 15])
hold on

% plot data & CI - pre, eyes open     
p1 = plot(x_closeup,  data_closeup_pre, 'Color', colours2(3, :), 'LineWidth', 2.5, 'LineStyle', ':');    
xlabel('frequency (Hz)')
ylabel('relative amplitude (µV)')
set(gca, 'FontSize', 16)
title('frontal region, eyes open', 'FontWeight', 'bold', 'FontSize', 18)
hold on
f1 = fill([x_closeup fliplr(x_closeup)],[data_closeup_pre + CI_closeup_pre fliplr(data_closeup_pre - CI_closeup_pre)], ...
    colours2(3, :), 'FaceAlpha', alpha, 'linestyle', 'none');
set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
hold on

% plot data & CI - post, eyes open    
p2 = plot(x_closeup,  data_closeup_post, 'Color', colours2(4, :), 'LineWidth', 2.5, 'LineStyle', ':');    
hold on
f2 = fill([x_closeup fliplr(x_closeup)],[data_closeup_post + CI_closeup_post fliplr(data_closeup_post - CI_closeup_post)], ...
    colours2(4, :), 'FaceAlpha', alpha, 'linestyle', 'none');
set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
hold on

% add legend 
legend([p1, p2], {'baseline' 'post medication'}, ...
    'FontSize', 16, 'location', 'southwest', 'edgecolor', [0.6, 0.6, 0.6])
hold off

% save figure
figure_name = 'rsEEG_alprazolam_closeup_1';
savefig(figure_name)
saveas(fig, [figure_name '.png'])
pause(3)

% update counter
figure_counter = figure_counter + 1;

% plot a close up - increase of eyes-closed low beta over the central region
% choose data
data_closeup_pre = data_pre_visual_closed(2, 8/header_high.xstep:14/header_high.xstep); data_closeup_post = data_post_visual_closed(2, 8/header_high.xstep:14/header_high.xstep);
CI_closeup_pre = CI_pre_visual_closed(2, 8/header_high.xstep:14/header_high.xstep); CI_closeup_post = CI_post_visual_closed(2, 8/header_high.xstep:14/header_high.xstep);

% adapt x
x_closeup = x_high(8/header_high.xstep:14/header_high.xstep);

% set limits of the figure
fig = figure(figure_counter);
plot(x_closeup, data_closeup_post + CI_closeup_post, 'b:', 'LineWidth', 0.5)
yl = get(gca, 'ylim');
clf
ylim([0, yl(2)])
xlim([12, 18])
hold on

% plot data & CI - pre, eyes closed     
p1 = plot(x_closeup,  data_closeup_pre, 'Color', colours2(3, :), 'LineWidth', 2.5);    
xlabel('frequency (Hz)')
ylabel('relative amplitude (µV)')
set(gca, 'FontSize', 16)
title('central region, eyes closed', 'FontWeight', 'bold', 'FontSize', 18)
hold on
f1 = fill([x_closeup fliplr(x_closeup)],[data_closeup_pre + CI_closeup_pre fliplr(data_closeup_pre - CI_closeup_pre)], ...
    colours2(3, :), 'FaceAlpha', alpha, 'linestyle', 'none');
set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
hold on

% plot data & CI - post, eyes closed  
p2 = plot(x_closeup,  data_closeup_post, 'Color', colours2(4, :), 'LineWidth', 2.5);    
hold on
f2 = fill([x_closeup fliplr(x_closeup)],[data_closeup_post + CI_closeup_post fliplr(data_closeup_post - CI_closeup_post)], ...
    colours2(4, :), 'FaceAlpha', alpha, 'linestyle', 'none');
set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
hold on

% add legend 
legend([p1, p2], {'baseline' 'post medication'}, ...
    'FontSize', 16, 'location', 'southwest', 'edgecolor', [0.6, 0.6, 0.6])
hold off

% save figure
figure_name = 'rsEEG_alprazolam_closeup_2';
savefig(figure_name)
saveas(fig, [figure_name '.png'])
pause(3)

% update counter
figure_counter = figure_counter + 1;

%% LOW FREQUENCIES: 1) calculate mean values
% preapre outcome variables
mean_low = [];
CI_low = [];
mean_bl_low = [];
CI_bl_low = [];

% mean values for each condition separately
for m = 1:datasize(1)
    for t = 1:datasize(2)
        for c = 1:datasize(3)
            for r = 1:datasize(5)
                for i = 1:datasize(6)
                    mean_low(m, t, c, r, i) = mean(data_low(m, t, c, :, r, i));
                    CI_low(m, t, c, r, i) = (std(data_low(m, t, c, :, r, i))/sqrt(datasize(4)))*z;
                end                
            end
        end
    end
end
datasize_mean = [size(mean_low, 1) size(mean_low, 2) size(mean_low, 3) size(mean_low, 4) size(mean_low, 5)];

% mean values for baseline conditions together 
for c = 1:datasize(3)
    for r = 1:datasize(5)
        for i = 1:datasize(6)
            data_bl = [squeeze(data_low(1, 1, c, :, r, i)); squeeze(data_low(2, 1, c, :, r, i))];
            mean_bl_low(c, r, i) = mean(data_bl);
            CI_bl_low(c, r, i) = (std(data_bl)/sqrt(length(data_bl)))*z;
        end                
    end
end

%% LOW FREQUENCIES: 2) plot baseline
% decide colour scheme
answer = questdlg('Would you like to set a new colour scheme?', 'Visualization - colours',...
    'Yes', 'No, load the preset', 'No, load the preset');
switch answer
    case 'Yes'
        for CS = 1:numel(ROI)
            colours(CS, :) = uisetcolor(['Choose colour for the ' ROI(CS).area ' ROI:']);
        end
    case 'No, load the preset'
        load('colours.mat');
end
clear answer

% determine x axis
x_start = header_low.xstart;
x_end = header_low.datasize(6) * header_low.xstep + header_low.xstart - header_low.xstep;
x_low = [x_start : header_low.xstep : x_end];

% ----- plot signal with closed eyes for all ROIs together -----
% choose data
data_visual_closed = squeeze(mean_bl_low(2, :, :));  

% set limits of the figure
fig = figure(figure_counter);
plot(x_low, data_visual_closed(1, :), 'b:', 'LineWidth', 0.5)
yl = get(gca, 'ylim');
clf
ylim([-0.2, yl(2)])
xlim([x_start, x_end])
hold on

% plot the data
for r = 1:size(data_visual_closed, 1)
    pl(r) = plot(x_low, data_visual_closed(r, :), 'Color', colours(r, :), 'LineWidth', 2);
    hold on
end
line([x_start, x_end], [0, 0], 'Color', [0, 0, 0], 'LineWidth', 0.5, 'LineStyle', ':');

% add name of frequency band
text((TOI(1).window(1) + TOI(1).window(2))/2, -0.09, TOI(1).sign, 'fontsize', 16, 'color', [0.6 0.6 0.6])
hold on

% add parameters  
xlabel('frequency (Hz)')
ylabel('relative amplitude (µV)')
set(gca, 'FontSize', 16)
title('eyes closed', 'FontWeight', 'bold', 'FontSize', 18)
hold on

% add legend 
legend(pl(1:5), {ROI(1:5).area}, 'FontSize', 16, 'position', [0.55, 0.55, 0.3, 0.3], 'edgecolor', [0.6 0.6 0.6])
hold off

% save figure
figure_name = ['rsEEG_baseline_closed_low'];
savefig(figure_name)
saveas(fig, [figure_name '.png'])
pause(3)

% update counter
figure_counter = figure_counter + 1;

% ----- plot signal with open eyes for all ROIs together -----
% choose data
data_visual_open = squeeze(mean_bl_low(1, :, :));  

% set limits of the figure
fig = figure(figure_counter);
plot(x_low, data_visual_closed(1, :), 'b:', 'LineWidth', 0.5)
yl = get(gca, 'ylim');
clf
ylim([-0.2, yl(2)])
xlim([x_start, x_end])
hold on

% plot the data
for r = 1:size(data_visual_closed, 1)
    pl(r) = plot(x_low, data_visual_closed(r, :), 'Color', colours(r, :), 'LineWidth', 2);
    hold on
end
line([x_start, x_end], [0, 0], 'Color', [0, 0, 0], 'LineWidth', 0.5, 'LineStyle', ':');

% add name of frequency band
text((TOI(1).window(1) + TOI(1).window(2))/2, -0.09, TOI(1).sign, 'fontsize', 16, 'color', [0.6 0.6 0.6])
hold on

% add parameters  
xlabel('frequency (Hz)')
ylabel('relative amplitude (µV)')
set(gca, 'FontSize', 16)
title('eyes open', 'FontWeight', 'bold', 'FontSize', 18)
hold on

% add legend 
legend(pl(1:5), {ROI(1:5).area}, 'FontSize', 16, 'position', [0.55, 0.55, 0.3, 0.3], 'edgecolor', [0.6 0.6 0.6])
hold off

% save figure
figure_name = ['rsEEG_baseline_open_low'];
savefig(figure_name)
saveas(fig, [figure_name '.png'])
pause(3)

% update counter
figure_counter = figure_counter + 1;

%% LOW FREQUENCIES: 3) plot difference pre X post
% decide colour scheme
answer = questdlg('Would you like to set a new colour scheme?', 'Visualization - colours',...
    'Yes', 'No, load the preset', 'No, load the preset');
switch answer
    case 'Yes'
        for CS = 1:4
            colours2(CS, :) = uisetcolor(['Choose colour ' num2str(CS)])
        end
    case 'No, load the preset'
        load('colours2.mat');
end
clear answer

% ----- plot pre X post medication for placebo -----
% choose data
data_pre_visual_open = squeeze(mean_low(1, 1, 1, :, :)); data_pre_visual_closed = squeeze(mean_low(1, 1, 2, :, :));  
CI_pre_visual_open = squeeze(CI_low(1, 1, 1, :, :)); CI_pre_visual_closed = squeeze(CI_low(1, 1, 2, :, :)); 
data_post_visual_open = squeeze(mean_low(1, 2, 1, :, :)); data_post_visual_closed = squeeze(mean_low(1, 2, 2, :, :));  
CI_post_visual_open = squeeze(CI_low(1, 2, 1, :, :)); CI_post_visual_closed = squeeze(CI_low(1, 2, 2, :, :)); 

% plot all regions
for r = 1:datasize(5)
    % set limits of the figure
    fig = figure(figure_counter);
    plot(x_low, data_post_visual_closed(r, :) + CI_post_visual_closed(r, :), 'b:', 'LineWidth', 0.5)
    yl = get(gca, 'ylim');
    clf
    ylim([-0.2, yl(2)])
    xlim([x_start, x_end])
    hold on
 
    % add names of frequency bands
    text((TOI(1).window(1) + TOI(1).window(2))/2, -0.09, TOI(1).sign, 'fontsize', 16, 'color', [0.6 0.6 0.6])
    
    % plot data & CI - pre, eyes open     
    p1 = plot(x_low,  data_pre_visual_open(r, :), 'Color', colours2(1, :), 'LineWidth', 2.5, 'LineStyle', ':');  
    line([x_start, x_end], [0, 0], 'Color', [0, 0, 0], 'LineWidth', 0.5, 'LineStyle', ':');
    xlabel('frequency (Hz)')
    ylabel('relative amplitude (µV)')
    set(gca, 'FontSize', 16)
    title(['placebo - ' ROI(r).area ' region'], 'FontWeight', 'bold', 'FontSize', 18)
    hold on
    f1 = fill([x_low fliplr(x_low)],[data_pre_visual_open(r, :) + CI_pre_visual_open(r, :) fliplr(data_pre_visual_open(r, :) - CI_pre_visual_open(r, :))], ...
        colours2(1, :), 'FaceAlpha', alpha, 'linestyle', 'none');
    set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
    hold on

    % plot data & CI - pre, eyes closed     
    p2 = plot(x_low,  data_pre_visual_closed(r, :), 'Color', colours2(1, :), 'LineWidth', 2.5);    
    hold on
    f2 = fill([x_low fliplr(x_low)],[data_pre_visual_closed(r, :) + CI_pre_visual_closed(r, :) fliplr(data_pre_visual_closed(r, :) - CI_pre_visual_closed(r, :))], ...
        colours2(1, :), 'FaceAlpha', alpha, 'linestyle', 'none');
    set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
    hold on
    
    % plot data & CI - post, eyes open   
    p3 = plot(x_low, data_post_visual_open(r, :), 'Color', colours2(2, :), 'LineWidth', 2.5, 'linestyle', ':');    
    hold on
    f3 = fill([x_low fliplr(x_low)],[data_post_visual_open(r, :) + CI_post_visual_open(r, :) fliplr(data_post_visual_open(r, :) - CI_post_visual_open(r, :))], ...
        colours2(2, :), 'FaceAlpha', alpha, 'linestyle', 'none');
    set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
    hold on
    
    % plot data & CI - alprazolam, eyes closed   
    p4 = plot(x_low,  data_post_visual_closed(r, :), 'Color', colours2(2, :), 'LineWidth', 2.5);    
    hold on
    f4 = fill([x_low fliplr(x_low)],[data_post_visual_closed(r, :) + CI_post_visual_closed(r, :) fliplr(data_post_visual_closed(r, :) - CI_post_visual_closed(r, :))], ...
        colours2(2, :), 'FaceAlpha', alpha, 'linestyle', 'none');
    set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
    hold on
    
    % add legend 
    legend([p1, p2, p3, p4], {'baseline - eyes open' 'baseline - eyes closed' 'post - eyes open' 'post - eyes closed'}, ...
        'FontSize', 16, 'position', [0.48, 0.55, 0.3 0.3], 'edgecolor', [0.6, 0.6, 0.6])
    hold off

    % save figure
    figure_name = ['rsEEG_placebo_' ROI(r).area '_low'];
    savefig(figure_name)
    saveas(fig, [figure_name '.png'])
    pause(3)

    % update counter
    figure_counter = figure_counter + 1;
end

% ----- plot pre X post medication for alprazolam -----
% choose data
data_pre_visual_open = squeeze(mean_low(2, 1, 1, :, :)); data_pre_visual_closed = squeeze(mean_low(2, 1, 2, :, :));  
CI_pre_visual_open = squeeze(CI_low(2, 1, 1, :, :)); CI_pre_visual_closed = squeeze(CI_low(2, 1, 2, :, :)); 
data_post_visual_open = squeeze(mean_low(2, 2, 1, :, :)); data_post_visual_closed = squeeze(mean_low(2, 2, 2, :, :));  
CI_post_visual_open = squeeze(CI_low(2, 2, 1, :, :)); CI_post_visual_closed = squeeze(CI_low(2, 2, 2, :, :)); 

% plot all regions
for r = 1:datasize(5)
    % set limits of the figure
    fig = figure(figure_counter);
    plot(x_low, data_post_visual_closed(r, :) + CI_post_visual_closed(r, :), 'b:', 'LineWidth', 0.5)
    yl = get(gca, 'ylim');
    clf
    ylim([-0.2, yl(2)])
    xlim([x_start, x_end])
    hold on
 
    % add names of frequency bands
    text((TOI(1).window(1) + TOI(1).window(2))/2, -0.09, TOI(1).sign, 'fontsize', 16, 'color', [0.6 0.6 0.6])
    
    % plot data & CI - pre, eyes open     
    p1 = plot(x_low,  data_pre_visual_open(r, :), 'Color', colours2(3, :), 'LineWidth', 2.5, 'LineStyle', ':');  
    line([x_start, x_end], [0, 0], 'Color', [0, 0, 0], 'LineWidth', 0.5, 'LineStyle', ':');
    xlabel('frequency (Hz)')
    ylabel('relative amplitude (µV)')
    set(gca, 'FontSize', 16)
    title(['alprazolam - ' ROI(r).area ' region'], 'FontWeight', 'bold', 'FontSize', 18)
    hold on
    f1 = fill([x_low fliplr(x_low)],[data_pre_visual_open(r, :) + CI_pre_visual_open(r, :) fliplr(data_pre_visual_open(r, :) - CI_pre_visual_open(r, :))], ...
        colours2(3, :), 'FaceAlpha', alpha, 'linestyle', 'none');
    set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
    hold on

    % plot data & CI - pre, eyes closed     
    p2 = plot(x_low,  data_pre_visual_closed(r, :), 'Color', colours2(3, :), 'LineWidth', 2.5);    
    hold on
    f2 = fill([x_low fliplr(x_low)],[data_pre_visual_closed(r, :) + CI_pre_visual_closed(r, :) fliplr(data_pre_visual_closed(r, :) - CI_pre_visual_closed(r, :))], ...
        colours2(3, :), 'FaceAlpha', alpha, 'linestyle', 'none');
    set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
    hold on
    
    % plot data & CI - post, eyes open   
    p3 = plot(x_low, data_post_visual_open(r, :), 'Color', colours2(4, :), 'LineWidth', 2.5, 'linestyle', ':');    
    hold on
    f3 = fill([x_low fliplr(x_low)],[data_post_visual_open(r, :) + CI_post_visual_open(r, :) fliplr(data_post_visual_open(r, :) - CI_post_visual_open(r, :))], ...
        colours2(4, :), 'FaceAlpha', alpha, 'linestyle', 'none');
    set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
    hold on
    
    % plot data & CI - alprazolam, eyes closed   
    p4 = plot(x_low,  data_post_visual_closed(r, :), 'Color', colours2(4, :), 'LineWidth', 2.5);    
    hold on
    f4 = fill([x_low fliplr(x_low)],[data_post_visual_closed(r, :) + CI_post_visual_closed(r, :) fliplr(data_post_visual_closed(r, :) - CI_post_visual_closed(r, :))], ...
        colours2(4, :), 'FaceAlpha', alpha, 'linestyle', 'none');
    set(get(get(f1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');  
    hold on
    
    % add legend 
    legend([p1, p2, p3, p4], {'baseline - eyes open' 'baseline - eyes closed' 'post - eyes open' 'post - eyes closed'}, ...
        'FontSize', 16, 'position', [0.48, 0.55, 0.3 0.3], 'edgecolor', [0.6, 0.6, 0.6])
    hold off

    % save figure
    figure_name = ['rsEEG_alprazolam_' ROI(r).area '_low'];
    savefig(figure_name)
    saveas(fig, [figure_name '.png'])
    pause(3)

    % update counter
    figure_counter = figure_counter + 1;
end


