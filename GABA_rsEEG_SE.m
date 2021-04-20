%% RS-EEG - COMPUTE SPECTRAL EXPONENT BETA 
% Written by Dominika for GABA-AD project (2021)
% Based on scripts and functions written by Michele A. Colombo from Massimini's lab, Milano
% 
% 1) Prepares data
%       - averages signals across ROIs
%       - cuts x and data according to target fband windows and saves in a
%         structure --> spect_exp.x; spect_exp.data
% 2) Fits the power using Michele's function fitPowerLaw3steps.m
%       - in 3 steps - first fit, alpha peak removal,second fit
%       - possible to visualize individual curves --> plot_i = 1
%       --> outcome variables:  intslo(1) = spect_exp.result.intercept
%                               intslo(2) = spect_exp.result.islope = SE beta
% 3) Performs group SE visualization
%       - extracts mean values 
%       - extracts individual SE change 
%       - visualize with box + scatter plots
% 4) Plots average log-log figures +- CI 
%       - calculates mean data, fits the second slope using fitPowerLaw3steps.m
%       - calculates CI = 95%, log-transforms
%       - plots on the logscale with peaks removed


%% parameters
clear all 
clc

% dataset
load('rsEEG_data_high.mat');
medication = {'placebo' 'alprazolam'}; 
time = {'pre' 'post'}; 
condition = {'open' 'closed'};
participant = 1:20;
xstep = 0.25; xoffset = 0.1;

% visualization
plot_i = 0;      
z = 1.96;
alpha = 0.2;
load('colours2.mat')
figure_counter = 1;

% outcome structure
spect_exp = struct;
spect_exp(1).band = 'low'; spect_exp(2).band = 'high'; spect_exp(3).band = 'broad'; 
spect_exp(1).window = [1 20]; spect_exp(2).window = [20 40]; spect_exp(3).window = [1 40]; 
spect_exp(1).method = 'ols'; spect_exp(2).method = 'ols'; spect_exp(3).method = 'ols';

% a = 2; m = 1; t = 1; c = 1; p = 1; r = 1; 

%% 1) prepare data
% average data across regions, square the amplitude --> power
for m = 1:size(data_high, 1)
    for t = 1:size(data_high, 2)
        for c = 1:size(data_high, 3)
            for p = 1:size(data_high, 4)
                for i = 1:size(data_high, 6)
                   data(m, t, c, p, i) = (mean(data_high(m, t, c, p, :, i)))^2;
                end
            end
        end
    end
end
clear data_high

% chop datasets by target fbands
for a = 1:numel(spect_exp)
    % crop data
    x_start = ceil((spect_exp(a).window(1) - xoffset) / xstep); 
    x_end = ceil((spect_exp(a).window(2) - xoffset) / xstep);
    spect_exp(a).data = data(:, :, :, :, x_start : x_end);   
    
    % indentify x
    x = spect_exp(a).window(1):xstep:spect_exp(a).window(2);
    spect_exp(a).x = x;
end
clear x_start x_end x data
clear a m t c p r i
   
%% 2) fit the power law 
% loop through target fband windows
for a = 1:numel(spect_exp)   
    % choose vector with frequency bins
    x = spect_exp(a).x;
    
    % loop through datasets
    for m = 1:size(spect_exp(a).data, 1)
        for t = 1:size(spect_exp(a).data, 2)
            for c = 1:size(spect_exp(a).data, 3)
                for p = 1:size(spect_exp(a).data, 4)
                    % launch figure if selected
                    if plot_i == 1
                        fig = figure(figure_counter)
                    end
                    
                    % call data
                    data = squeeze(spect_exp(a).data(m, t, c, p, :))';
                    
                    % choose plot colour
                    col = colours2((m - 1)*2 + t, :);
                    
                    % fit current dataset
                    [intslo, stats, amps] = fitPowerLaw3steps(x, data, spect_exp(a).method, plot_i, col); 
                    
                    % fill in the outcome structure
                    spect_exp(a).result.intercept(m, t, c, p) = intslo(1);
                    spect_exp(a).result.slope(m, t, c, p) = intslo(2);  
                    
                    % update counter if selected
                    if plot_i == 1
                        figure_counter = figure_counter + 1;
                    end
                end 
            end
        end
    end
end
clear dataset_name x data col intslo stats amps fig 
clear a m t c p i

% save outcome structure
save('rsEEG_spect_exp.mat', 'spect_exp')

%% 3) group visualization 
% ----- extract mean SE values -----
for a = 1:numel(spect_exp)  
    for m = 1:size(spect_exp(a).data, 1)
        for t = 1:size(spect_exp(a).data, 2)
            for c = 1:size(spect_exp(a).data, 3)
                avg_se(m, t, c) = mean(spect_exp(a).result.slope(m, t, c, :));                
                avg_sem(m, t, c) = std(spect_exp(a).result.slope(m, t, c, :))/sqrt(length(participant));                               
            end
        end
    end
end
clear a m t c

% ----- plot group boxplot -----
for a = 1:numel(spect_exp)  
    for c = 1:length(condition)
        % choose the data
        data_visual = [];
        for m = 1:length(medication)
            for t = 1:length(time)
                data_i = squeeze(spect_exp(a).result.slope(m, t, c, :));
                data_visual = cat(2, data_visual, data_i);
            end
        end

        % plot group boxplot
        fig = figure(figure_counter);        
        boxplot(data_visual, 'color', colours2)
        hold on
        
        % plot the lines - placebo
        for p = 1:length(participant)
            p_placebo(p) = plot([1 2], data_visual(p, [1 2]), '-o',...
                'Color', [0.75, 0.75, 0.75],...
                'MarkerSize', 10,...
                'MArkerEdge', 'none');
            hold on
        end
        
        % plot the lines - alprazolam
        for p = 1:length(participant)
            p_alprazolam(p) = plot([3 4], data_visual(p, [3 4]), '-',...
                'Color', [0.75, 0.75, 0.75],...
                'MarkerSize', 10,...
                'MArkerEdge', 'none');
            hold on
        end

        % plot the markers
        for b = 1:4
            scat(b) = scatter(repelem(b, length(participant)), data_visual(:, b),...
                75, colours2(b, :), 'filled');
            hold on
        end
        
        % add parameters
        figure_title = ['SPECTRAL EXPONENT: ' spect_exp(a).band ' band window, eyes ' condition{c}];
        figure_name = ['rsEEG_' spect_exp(a).band '_' condition{c}];
        yl = get(gca, 'ylim');
        set(gca, 'xtick', 1:4, 'xticklabel', {'pre' 'post' 'pre' 'post'})
        set(gca, 'Fontsize', 16)
        title(figure_title, 'FontWeight', 'bold', 'FontSize', 18)
        xlabel('time relative to medication'); ylabel('spectral exponent');
        ylim([yl(1) - 0.05, yl(2) + 0.1])
        hold on
        
        % add text
        for m = 1:length(medication)
            for t = 1: length(time)
                txt((m - 1)*2 + t) = text((m - 1)*2 + t - 0.3, yl(1) - 0.01,...
                    ['mean: ' num2str(round(avg_se(m, t, c), 3))],...
                    'fontsize', 12, 'color', [0.4 0.4 0.4]);
                hold on
            end
        end
        txt(end + 1) = text(1.3, yl(2) + 0.05, 'placebo', 'fontsize', 16, 'color', colours2(2, :));
        txt(end + 1) = text(3.3, yl(2) + 0.05, 'alprazolam', 'fontsize', 16, 'color', colours2(4, :));
        hold off

        % save the figure
        pause(5)        
        savefig([figure_name '.fig'])
        saveas(fig, [figure_name '.png'])

        % update the counter
        figure_counter = figure_counter + 1;        
    end
end
clear data_i data_visual figure_name figure_title p_placebo p_alprazolam scat yl txt 
clear a b m c p

% ----- calculate individual change ----- 
for a = 1:numel(spect_exp)  
    for m = 1:size(spect_exp(a).data, 1)
        for c = 1:size(spect_exp(a).data, 3)
            for p = 1:size(spect_exp(a).data, 4)
                spect_exp(a).result.slopechange(m, c, p) = spect_exp(a).result.slope(m, 2, c, p) - spect_exp(a).result.slope(m, 1, c, p);                                             
            end
        end
    end
end
save('rsEEG_spect_exp.mat', 'spect_exp')
clear a m c p

% ----- plot group boxplot of individual change ----- 
for a = 1:numel(spect_exp)  
    % choose the data
    data_visual = [];
    for c = 1:length(condition) 
        for m = 1:length(medication)
            data_i = squeeze(spect_exp(a).result.slopechange(m, c, :));
            data_visual = cat(2, data_visual, data_i);
        end
    end
     
    % plot group boxplot
    col = colours2([2 4 2 4], :);
    fig = figure(figure_counter);        
    boxplot(data_visual, 'color', col)
    hold on
    
    % add zero line
    xl = get(gca, 'xlim');
    h = line([xl(1) xl(2)], [0, 0], 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 2);
    hold on

    % plot the lines - eyes open
    for p = 1:length(participant)
        p_open(p) = plot([1 2], data_visual(p, [1 2]), '-o',...
            'Color', [0.75, 0.75, 0.75],...
            'MarkerSize', 10,...
            'MArkerEdge', 'none');
        hold on
    end

    % plot the lines - eyes closed
    for p = 1:length(participant)
        p_closed(p) = plot([3 4], data_visual(p, [3 4]), '-',...
            'Color', [0.75, 0.75, 0.75],...
            'MarkerSize', 10,...
            'MArkerEdge', 'none');
        hold on
    end

    % plot the markers
    for b = 1:4
        scat(b) = scatter(repelem(b, length(participant)), data_visual(:, b),...
            75, col(b, :), 'filled');
        hold on
    end

    % add parameters
    figure_title = ['SPECTRAL EXPONENT CHANGE - ' spect_exp(a).band ' band window'];
    figure_name = ['rsEEG_SEchange_' spect_exp(a).band];
    yl = get(gca, 'ylim');
    set(gca, 'xtick', 1:4, 'xticklabel', {'placebo' 'alprazolam' 'placebo' 'alprazolam'})
    set(gca, 'Fontsize', 16)
    title(figure_title, 'FontWeight', 'bold', 'FontSize', 18)
    xlabel('medication'); ylabel('spectral exponent change');
    ylim([yl(1) - 0.05, yl(2) + 0.1])
    hold on

    % add text
    for b = 1:4
        txt(b) = text(b - 0.3, yl(1) - 0.01,...
            ['mean: ' num2str(round(mean(data_visual(:, b)), 3))],...
            'fontsize', 12, 'color', [0.4 0.4 0.4]);
        hold on
    end
    txt(end + 1) = text(1.25, yl(2) + 0.05, 'eyes open', 'fontsize', 16, 'color', [0 0 0]);
    txt(end + 1) = text(3.25, yl(2) + 0.05, 'eyes closed', 'fontsize', 16, 'color', [0 0 0]);
    hold off

    % save the figure
    pause(5)        
    savefig([figure_name '.fig'])
    saveas(fig, [figure_name '.png'])

    % update the counter
    figure_counter = figure_counter + 1;        

end
clear data_i data_visual xl h figure_name figure_title p_open p_closed scat yl txt col fig
clear a b m c p

%% 4) plot average log-log figures
% separate figure per fband window, medication, and condition
for a = 1:numel(spect_exp)  
    % choose vector with frequency bins
    x = spect_exp(a).x;
    
    % loop through datasets
    for c = 1:length(condition) 
        % launch a figure
        fig = figure(figure_counter);
        
        % add both med and time conditions
        for  m = 1:length(medication)
            
            subplot(1, 2, m)
            for t = 1:length(time)
                % ----- calculate mean values -----
                for i = 1:size(spect_exp(a).data, 5)
                    spect_exp(a).avgdata.data(m, t, c, i) = mean(squeeze(spect_exp(a).data(m, t, c, :, i)));
                    spect_exp(a).avgdata.CI(m, t, c, i, 1) = spect_exp(a).avgdata.data(m, t, c, i) + (std(squeeze(spect_exp(a).data(m, t, c, :, i)))/sqrt(length(participant)))*z;     
                    spect_exp(a).avgdata.CI(m, t, c, i, 2) = spect_exp(a).avgdata.data(m, t, c, i) - (std(squeeze(spect_exp(a).data(m, t, c, :, i)))/sqrt(length(participant)))*z; 
                end
                
                 % log-transform CI
                 for i = [1 2]
                     y = squeeze(spect_exp(a).avgdata.CI(m, t, c, :, i))';
                     [x_int, spect_exp(a).avgdata.CI_int(m, t, c, :, i)] = log_int(x, y);
                 end

                % ----- fit power law to mean values -----
                % choose data
                data = squeeze(spect_exp(a).avgdata.data(m, t, c, :))';
                
                % choose plot colour
                col = colours2((m - 1)*2 + t, :);

                % fit current dataset
                [intslo, stats, amps, devs] = fitPowerLaw3steps(x, data, spect_exp(a).method, 0, col); 
                
                % fill in the outcome structure
                spect_exp(a).avgdata.data_int(m, t, c, :) = amps.obs;
                spect_exp(a).avgdata.x_int(m, t, c, :) = amps.frex;   
                spect_exp(a).avgdata.intercept(m, t, c) = intslo(1);
                spect_exp(a).avgdata.slope(m, t, c) = intslo(2);  
                
                % ----- plot the final signal -----                
                % plot the original signal
                pl_dot((m - 1)*2 + t) = loglog(amps.frex, amps.obs, ':','color', col, 'linewidth', 1.5); 
                hold on
                
                % highlight kept points
                x_kept = amps.frex; 
                x_kept(devs.rej)= nan;
                pl_full((m - 1)*2 + t) = loglog(x_kept, amps.obs, '-', 'color', col, 'LineWidth', 2); 
                hold on;
                
                % add CI shading
                ci_visual = squeeze(spect_exp(a).avgdata.CI_int(m, t, c, :, :))';
                shade((m - 1)*2 + t) = fill([amps.frex fliplr(amps.frex)],[ci_visual(1, :) fliplr(ci_visual(2, :))], ...
                    col, 'FaceAlpha', alpha, 'linestyle', 'none');
                hold on

                % plot second fitted line 
                lfit((m - 1)*2 + t) = plot(amps.frex, amps.pred, '--', 'color', [0.4 0.4 0.4], 'LineWidth', 1.5);
                hold on
            end
            % add other parameters
            xlim([amps.frex(1) amps.frex(end)])
            ylim([0, amps.pred(1) + 0.2])
            xlabel('Hz'); ylabel('relative amplitude');
            set(gca, 'FontSize', 16)            
            
            % add legend and text
            sgtitle(['\fontsize{18}' spect_exp(a).band ' frequencies - eyes ' condition{c}])
            title(['\fontsize{16}' medication{m}])
            legend([pl_full((m-1)*2 + 1) pl_full((m-1)*2 + 2)], {'baseline', 'post medication'},...
                'FontSize', 16, 'location', 'northeast', 'edgecolor', [0.8 0.8 0.8])
            hold on
                      
        end        
        % save figure  
        pause(10)    
        figure_name = ['rsEEG_SE_' spect_exp(a).band '_' condition{c}];
        savefig([figure_name '.fig'])
        saveas(fig, [figure_name '.png'])

        % update the counter
        figure_counter = figure_counter + 1;             
    end
end
clear x y x_int data col x_kept intslo stats amps devs ci_visual pl_dot pl_full shade lfit fig figure_name
clear a m t c i















