%% RS-EEG - EXTRACTION OF MEAN FREQUENCY BAND AMPLITUDE 
% Written by Dominika for GABA-AD project (2021)
% 
% ----- AMPLITUDES -----
% 1) Prepares data
%       - Loads data, pools eletrodes in areas according to predefined ROIs
%       - Merges data of all participants into a 6D -mat matrix, saves 
% 2) Identifies individual alpha frequency (IAF) 
%       - In eyes-closed datasets
%       - Allows user to visualize signal from all ROIs and manually select
%         individual time window for IAF search --> saves as 'IAF_window.mat'
%       - Identifies individual IAF as the frequency within the corresponding 
%         time window with the highest amplitude 
%       --> saves as table - 'IAF.mat'; 'IAF.csv'   
% 3) Plots group IAF values
%       - Group mean values - boxplot --> saves as 'IAF_<ROI>.fig + png'
%       - Individual IAF change - paired scatterplot
%       --> saves values as 'IAF_change.mat' and figuers as 'IAF_scatter_<ROI>.fig + png'
% 4) Identifies individual transition frequency (TF)
%       - In baseline datasets - eyes-closed vs eyes-open
%       - Allows user to visualize signal from occipital ROI and manually 
%         select TF (separately for each medication) 
%       --> saves as table - 'TF.mat'; 'TF.csv'   
% 5) Calculates individual limits of frequency bands based on IAF and TF
%       - rules predefined in TOI (struct) - see parameters
%       --> saves as table - 'fband.mat'; 'fband.csv' 
%       - Plot fband distribution based on IAF and TF in data from 
%         an example participant
% 6) Extracts mean amplitude values 
%       - separately from all frequency bands (=TOIs) and all ROIs
%       --> saves as table - 'fband_amplitude.mat'; 'fband_amplitude.csv'
% 7) Plots mean amplitudes 
%       - barplots for each ROI
%       - normalizes amplitudes to baseline = change in %
%       - line plot for each frequency subband --> visual identification of
%       effect of the medication
% 8) Extracts beta band amplitude
%       - plots changes in beta mean amplitude
%       - ROI: central region, TOI: beta1, beta2, broad band beta (beta1 + beta2)
%       --> saves values as 'beta.mat' 
% 9) Calculates alpha attenuation coefficient (AAC)
%       - for individual alpha subbands + broad alpha band
%       - calculates individual AAC and AAC change 
%       - plots box + scatter plots, saves figures
%       --> saves values as 'AAC.mat' and 'AAC.csv'
% 
% ----- SPECTRAL EXPONENT -----
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

%% AMPLITUDES: parameters
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
TOI(1).band = 'delta'; TOI(1).definition = 'window(1) = 0.1; window(2) = tf - 2;'; TOI(1).sign = '\delta'; 
TOI(2).band = 'theta'; TOI(2).definition = 'window(1) = tf - 2; window(2) = tf;'; TOI(2).sign = '\theta'; 
TOI(3).band = 'alpha1'; TOI(3).definition = 'window(1) = tf; window(2) = tf + (iaf - tf)/2;'; TOI(3).sign = '\alpha 1'; 
TOI(4).band = 'alpha2'; TOI(4).definition = 'window(1) =  tf + (iaf - tf)/2; window(2) = iaf;'; TOI(4).sign = '\alpha 2'; 
TOI(5).band = 'alpha3'; TOI(5).definition = 'window(1) = iaf; window(2) = iaf + (iaf - tf)/2;'; TOI(5).sign = '\alpha 3'; 
TOI(6).band = 'beta1'; TOI(6).definition = 'window(1) = iaf + (iaf - tf)/2; window(2) = 20;'; TOI(6).sign = '\beta 1'; 
TOI(7).band = 'beta2'; TOI(7).definition = 'window(1) = 20; window(2) = 30;'; TOI(7).sign = '\beta 2'; 
TOI(8).band = 'gamma'; TOI(8).definition = 'window(1) = 30; window(2) = 45;'; TOI(8).sign = '\gamma'; 

% import header basis
load([prefix_1 ' high ' prefix_2 '1 placebo pre EEGcont open.lw6'], '-mat')
header_high = header;
load([prefix_1 ' low ' prefix_2 '1 placebo pre EEGcont open.lw6'], '-mat')
header_low = header;
clear header

% statistics
z = 1.96;

% visualization
figure_counter = 1;
col = [1, 0.55, 0.55];
load('colours.mat');

% in case of repeted use:
load('rsEEG_data_high.mat'); load('rsEEG_data_low.mat') 
load('IAF.mat'); load('TF.mat')
load('fband.mat'); load('rsEEG_fband_amplitude.mat')

%% 1) load data, pool channels into ROIs

% ----- HIGHER FREQUENCIES -----
% load individual subject data and calculate mean signal for each ROI
data_high = [];
for m = 1:numel(medication)
    for t = 1:numel(time)
        for c = 1:numel(condition)
            for p = 1:numel(participant)
                % load the dataset
                dataset_name = [prefix_1 ' high ' prefix_2 num2str(participant(p)) ' ' medication{m} ' ' time{t} ' EEGcont ' condition{c} '.mat'];
                load(dataset_name)
                
                % calculate mean signal for each ROI area                
                for r = 1:numel(ROI)                    
                    % load signal from target electrodes
                    channels2merge = [];
                    for e = 1:numel(ROI(r).electrodes)
                        % identify target electrodes
                        position = find(contains(labels, ROI(r).electrodes{e}) & cellfun('length', labels) == length(ROI(r).electrodes{e}));
                        channels2merge(e, :) = squeeze(data(:, position, :, :, :, :))';                        
                    end
                    
                    % calculate mean signal
                    for i = 1:length(channels2merge)
                        data_high(m, t, c, p, r, i) = mean(channels2merge(:, i));
                    end
                end
            end
        end
    end
end

% verify the datasize
datasize_high = [size(data_high, 1) size(data_high, 2) size(data_high, 3) size(data_high, 4) size(data_high, 5) size(data_high, 6)];
disp('High frequencies: data matrix successfully created. Datsize:')
disp(datasize_high)

% ----- LOWER FREQUENCIES -----
% load individual subject data and calculate mean signal for each ROI
data_low = [];
for m = 1:numel(medication)
    for t = 1:numel(time)
        for c = 1:numel(condition)
            for p = 1:numel(participant)
                % load the dataset
                dataset_name = [prefix_1 ' low ' prefix_2 num2str(participant(p)) ' ' medication{m} ' ' time{t} ' EEGcont ' condition{c} '.mat'];
                load(dataset_name)
                
                % calculate mean signal for each ROI area                
                for r = 1:numel(ROI)                    
                    % load signal from target electrodes
                    channels2merge = [];
                    for e = 1:numel(ROI(r).electrodes)
                        % identify target electrodes
                        position = find(contains(labels, ROI(r).electrodes{e}) & cellfun('length', labels) == length(ROI(r).electrodes{e}));
                        channels2merge(e, :) = squeeze(data(:, position, :, :, :, :))';                        
                    end
                    
                    % calculate mean signal
                    for i = 1:length(channels2merge)
                        data_low(m, t, c, p, r, i) = mean(channels2merge(:, i));
                    end
                end
            end
        end
    end
end

% verify the datasize
datasize_low = [size(data_low, 1) size(data_low, 2) size(data_low, 3) size(data_low, 4) size(data_low, 5) size(data_low, 6)];
disp('Low frequencies: data matrix successfully created. Datsize:')
disp(datasize_low)

clear dataset_name position channels2merge data

% % save data for future use
% save('rsEEG_data_high.mat', 'data_high');
% save('rsEEG_data_low.mat', 'data_low');


%% 2) identify individual IAF

% ----- CHOOSE LOOK-UP WINDOW -----
% parameters of visualization
window = [6, 14];
x = [window(1) : header_high.xstep : window(2)];
x_start = ceil((window(1)-header_high.xstart)/header_high.xstep);
x_end = ceil((window(2)-header_high.xstart)/header_high.xstep);
def_window_start = 9; 
def_window_span = 2; 

% loop through subjects
window_final = [];
for p = 1:length(participant)
    % choose data for visualization 
    data_visual = cat(1, squeeze(data_high(1, 1, 2, p, :, x_start:x_end)), squeeze(data_high(2, 1, 2, p, :, x_start:x_end)));

    % check if vector size matches
    if length(data_visual(1,:)) ~= length(x)
        diff = length(data_visual(1,:)) - length(x);
        if diff > 0
            data_visual = data_visual(:, 1:end - diff);
        elseif diff < 0
            data_visual = data_visual(:, 1:end + diff);
        end
    end
    
    % set default parameters
    window_span = def_window_span;
    window_start = def_window_start;
    
    % identify the final window for current participant
    finish = 0;
    while finish == 0;
        % launch the figure
        fig = figure(figure_counter);
        
        % plot the default window 
        plot(x, data_visual(:, :), 'k:', 'LineWidth', 0.5)
        yl = ylim;
        clf
        xlim([window(1), window(2)])
        ylim(yl)
        hr = rectangle('Position', [window_start, yl(1), window_span, yl(2)-yl(1)], 'FaceColor', col, 'EdgeColor', 'none');
        title(['IAF window - participant ' num2str(participant(p))], 'FontSize', 16)
        set(gca, 'FontSize', 16)
        hold on
        
        % plot the data timeseries
        plot(x, data_visual(:, :), 'k')
        hold on
        
        % choose the peak position and plot it
        pos_x = get_position(gca);  
        hl = line([pos_x, pos_x], [yl(1), yl(2)], 'Color', 'red', 'LineWidth', 2); 
        
        % replot the window
        window_start = floor(pos_x/header_high.xstep)*header_high.xstep - window_span/2;
        set(hr, 'Position', [window_start, yl(1), window_span, yl(2)-yl(1)]);
        
        % check if OK to leave
        answer = questdlg('Do you want to proceed?', ['IAF window - participant ' num2str(participant(p))],...
            'Yes, save window values.', 'No, I want to adjust center.', 'No, I want to adjust span.', 'Yes, save window values.');
        switch answer
            case 'Yes, save window values.'
                % close the figure
                close(fig)

                % exit the while loop
                finish = 1;

            case 'No, I want to adjust center.'
            case 'No, I want to adjust span.'
                % set the new span value
                prompt = {'New span in x units:'};
                dlgtitle = 'Set the window span';
                dims = [1 35];
                definput = {'2'};
                window_span = str2num(cell2mat(inputdlg(prompt,dlgtitle,dims,definput)));
                clear answer2 prompt dlgtitle dims definput                      
        end
    end

    % save window parameters
    window_final(p, 1) = window_start; 
    window_final(p, 2) = window_span; 
end
save('IAF_window.mat', 'window_final') 
clear answer finish window window_start window_span final_span def_window_span def_window_start data_visual hl hr pos_x x x_start x_end 


% ----- IDENTIFY IAF (EYES CLOSED) -----
% create outcome table
IAF = table;
IAF.subject = (repelem(participant, length(medication) * length(time)))';
IAF.medication = repmat((repelem(medication, length(time)))', length(participant), 1);
IAF.time = repmat(time', length(participant) * length(medication), 1);

% loop through ROIs
for r = 1:5
    % create outcome column
    statement = [ROI(r).area ' = [];'];
    eval(statement)
        
    % loop through datasets
    for p = 1:length(participant)
        for m = 1:length(medication)
            for t = 1:length(time)  
                % choose data
                x_start = ceil((window_final(p, 1) - header_high.xstart)/header_high.xstep);
                x_end = ceil(((window_final(p, 1) + window_final(p, 2)) - header_high.xstart)/header_high.xstep);
                data_i = squeeze(data_high(m, t, 2, p, r, x_start:x_end));
                
                % calculate iaf
                index = x_start + find(data_i == max(data_i)) - 1;
                iaf_i = index * header_high.xstep;                          % no need to add xstart, 0.1 is under the resolution 
                
                % append the result
                statement = [ROI(r).area ' = [' ROI(r).area '; iaf_i];'];
                eval(statement)                
            end
        end
    end
    statement = ['IAF.' ROI(r).area ' = ' ROI(r).area ';'];
    eval(statement)
end

% save outcome table
save('IAF.mat', 'IAF')
writetable(IAF, 'IAF.csv')

clear statement data_i x_start x_end index iaf_i


%% 3) visualize mean IAF 

% ----- IAF MEAN VALUES -----        
% extract mean values
for m = 1:length(medication)
    for t = 1:length(time)
        for r = 1:numel(ROI)
            % choose data
            rows = (categorical(IAF.medication) == medication{m} & categorical(IAF.time) == time{t});
            data_i = table2array(IAF(rows, 3 + r));
            
            % calculate mean and CI (95%)
            IAF_mean(m, t, r) = mean(data_i);
            IAF_CI(m, t, r) = (std(IAF_mean(m, t, r))/sqrt(length(participant)))*z;
        end
    end 
end
clear rows data_i

% plot a box plot for each ROI
for r = 1:numel(ROI)
    % choose the data
    data_visual = [];
    for m = 1:length(medication)
        for t = 1:length(time)
            rows = (categorical(IAF.medication) == medication{m} & categorical(IAF.time) == time{t});
            data_visual = cat(2, data_visual, table2array(IAF(rows, 3 + r)));
        end
    end
        
    % launch the figure
    fig = figure(figure_counter);
    boxplot(data_visual, 'colors', colours2([2 2 4 4], :))
    set(gca, 'xtick', 1:4, 'xticklabel', {'pre' 'post' 'pre' 'post'})
    set(gca, 'Fontsize', 16)
    title([ROI(r).area ' region'], 'FontWeight', 'bold', 'FontSize', 18)
    
    % save the figure
    figure_name = ['IAF_' ROI(r).area];
    savefig([figure_name '.fig'])
    saveas(fig, [figure_name '.png'])
    
    % update the counter
    figure_counter = figure_counter + 1;
end
clear rows data_visual figure_name

% plot box plots for baseline IAF
% choose the data
rows = (categorical(IAF.time) == time{1});
data_visual = table2array(IAF(rows, 4:end));

% launch the figure
fig = figure(figure_counter);
boxplot(data_visual, 'colors', colours)
set(gca, 'xtick', 1:5, 'xticklabel', {ROI(:).area})
set(gca, 'Fontsize', 16)
title('PAF across topographical regions', 'FontWeight', 'bold', 'FontSize', 18)
xlabel('region')
ylabel('mean PAF (Hz)')

% save the figure
figure_name = 'IAF_baseline';
savefig([figure_name '.fig'])
saveas(fig, [figure_name '.png'])

% update the counter
figure_counter = figure_counter + 1;
clear rows data_visual figure_name


% ----- IAF CHANGE -----
% extract individual IAF change
for p = 1:length(participant)
    for m = 1:length(medication)
        for r = 1:numel(ROI)  
            % choose data
            rows = (IAF.subject == participant(p) & categorical(IAF.medication) == medication{m});
            data_i = table2array(IAF(rows, 3 + r));
            
            % calculate IAF change
            IAF_change(p, m, r) = data_i(2) - data_i(1);
        end 
    end
end
save('IAF_change.mat', 'IAF_change')
clear rows data_i    

% plot a line graph for each region
for r = 1:numel(ROI)
    fig = figure(figure_counter);
    % plot the lines
    for p = 1:length(participant)
        pl(p) = plot([1 2], squeeze(IAF_change(p, :, r)), '-ko',...
            'MarkerSize', 10,...
            'MArkerEdge', 'none');
        hold on
    end
    
    % add zero line
    h = line([0.75 2.25], [0, 0], 'Color', [0.75, 0.75, 0.75], 'LineStyle', ':', 'LineWidth', 2);
    hold on
    
    % plot the markers
    scattercolor = colours2([2 4], :);
    for m = 1:length(medication)
        sc(m) = scatter(repelem(m, length(participant)), squeeze(IAF_change(:, m, r)),...
            75, scattercolor(m, :), 'filled');
        hold on
    end
       
    % adjust parameters
    xlim([0.75 2.25])
    ylim([-3 3])
    set(gca, 'xtick', 1:2, 'xticklabel', medication)
    set(gca, 'Fontsize', 14)
    title([ROI(r).area ' region'], 'FontWeight', 'bold', 'FontSize', 16)
    xlabel('medication')
    ylabel('IAF change (Hz)')
    
    % save the figure
    figure_name = ['IAF_scatter_' ROI(r).area];
    savefig([figure_name '.fig'])
    saveas(fig, [figure_name '.png'])
    
    % update figure counter
    figure_counter = figure_counter + 1;
end
clear pl sc

%% 4) identify individual TF
% create outcome variables
TF = table;
TF.subject = (repelem(participant, length(medication)))';
TF.medication = repmat((medication)', length(participant), 1);
tf_occipital = []; 

% parameters of visualization
window = [3, 15];
x = [window(1) : header_high.xstep : window(2)];
x_start = ceil((window(1)-header_high.xstart)/header_high.xstep);
x_end = ceil((window(2)-header_high.xstart)/header_high.xstep);

% loop through subjects - data from occipital region
for p = 1:length(participant)
    % launch the figure
    fig = figure(figure_counter);
    
    % loop through sessions
    for m = 1:length(medication)
        % choose baseline data for visualization 
        data_visual = cat(1, squeeze(data_high(m, 1, 1, p, 5, x_start:x_end))', squeeze(data_high(m, 1, 2, p, 5, x_start:x_end))');

        % check if vector size matches
        if length(data_visual(1,:)) ~= length(x)
            diff = length(data_visual(1,:)) - length(x);
            if diff > 0
                data_visual = data_visual(:, 1:end - diff);
            elseif diff < 0
                data_visual = data_visual(:, 1:end + diff);
            end
        end

        % identify TF manually
        finish = 0;
        while finish == 0;
            % plot the datasets for current medication
            subplot(2, 4, [(m - 1)*2 + 1, (m - 1)*2 + 2, (m - 1)*2 + 5, (m - 1)*2 + 6])
            for c = 1:size(data_visual, 1)
                if c == 1
                    pl(c) = plot(x, data_visual(c, :), 'color', colours2((m - 1)*2 + 1, :), 'linewidth', 2, 'linestyle', ':');
                elseif c == 2
                    pl(c) = plot(x, data_visual(c, :), 'color', colours2((m - 1)*2 + 1, :), 'linewidth', 2);
                end
                hold on
            end
            
            % add parameters
            yl = get(gca,'ylim');
            xlim([window(1), window(2)])
            ylim(yl)
            title(medication{m}, 'FontSize', 16)
            set(gca, 'FontSize', 14)
            xlabel('frequency (Hz)')
            ylabel('relative amplitude (�V)')
            hold on
            
            % add IAF
            row = (IAF.subject == participant(p) & categorical(IAF.medication) == medication{m} & categorical(IAF.time) == time{1});
            x_iaf = IAF.occipital(row);
            hl_iaf = line([x_iaf, x_iaf], [yl(1), yl(2)], 'color', [0, 0, 0], 'LineWidth', 1.5); 
            hold on

            % choose TF and plot it
            x_tf = get_position(gca);  
            hl_tf = line([x_tf, x_tf], [yl(1), yl(2)], 'Color', [0, 0, 0], 'LineWidth', 1.5, 'LineStyle', '--'); 
            hold on

           % check if OK to leave
            answer = questdlg('Do you want to proceed?', ['TF - participant ' num2str(participant(p)) ', ' medication{m}],...
                'Yes, save TF.', 'No, I want to make adjustments.', 'Yes, save TF.');
            switch answer
                case 'Yes, save TF.'
                    % exit the while loop
                    finish = 1;

                case 'No, I want to make adjustments.'                 
            end
        end
        
        % calculate final TF value
        x_tf = ceil(x_tf/header_high.xstep) * header_high.xstep;

        % append TF
        tf_occipital = [tf_occipital; x_tf]; 
    end
    
    % add main title
    suptitle(['participant ' num2str(participant(p)) ' - occipital region'])
    
    % save figure
    figure_name = ['TF_YC' num2str(participant(p))];
    savefig([figure_name '.fig'])
    saveas(fig, [figure_name '.png'])
    
    %update the counter
    figure_counter = figure_counter +1;
end

% play a celebratory sound at the end of the task <3
tune = load('handel.mat');
sound(tune.y, tune.Fs)

% complete the outcome table and save
TF.occipital = tf_occipital;
save('TF.mat', 'TF') 
writetable(TF, 'TF.csv')

clear window x x_start x_end answer finish window data_visual row x_iaf hl_iaf x_tf hl_tf tf_occipital tune

%% 5) calculate individual limits of frequency bands 

% ----- EXTRACT DATA -----
% based on occipital IAF and TF
% create outcome table
fband = table;
fband.subject = (repelem(participant, length(medication)))';
fband.medication = repmat(medication', length(participant), 1);

% loop through frequency bands
for f = 1:numel(TOI)
    % create outcome column
    statement = [TOI(f).band ' = [];'];
    eval(statement) 
    
    % loop through participants and session
    for p = 1:length(participant)
        for m = 1:length(medication)
            % call iaf 
            row = (IAF.subject == participant(p) & categorical(IAF.medication) == medication{m} & categorical(IAF.time) == time{1});
            iaf = IAF.occipital(row);
            
            % call tf
            row = (TF.subject == participant(p) & categorical(TF.medication) == medication{m});
            tf = TF.occipital(row);
                
            % calculate fband boundaries
            statement = TOI(f).definition;
            eval(statement) 

            % append the result
            statement = [TOI(f).band ' = [' TOI(f).band '; {window}];'];
            eval(statement)              
        end
    end 

    % complete the outcome table and save
    statement = ['fband.' TOI(f).band ' = ' TOI(f).band ';'];   
    eval(statement)   
end

% save the outcome table
save('fband.mat', 'fband')
writetable(fband, 'fband.csv')

clear statement row iaf tf window delta theta alpha1 alpha2 alpha3 beta1 beta2 gamma tune

%  ----- PLOT REPRESENTATIVE PARTICIPANT -----
% choose data
data_visual_open = squeeze(data_high(1, 1, 1, 1, 5, [8 : end-1]))'; 
data_visual_closed = squeeze(data_high(1, 1, 2, 1, 5, [8 : end-1]))';  

% x axis
window_visual = [2 45];
x = [window_visual(1) : header_high.xstep : window_visual(2)];

% logtransform 
[x_log, open_log] = log_int(x, data_visual_open);
[x_log, closed_log] = log_int(x, data_visual_closed);

% identify IAF and TF
rows = (IAF.subject == 1 & categorical(IAF.medication) == 'placebo' & categorical(IAF.time) == 'pre');
iaf = IAF.occipital(rows);
rows = (TF.subject == 1 & categorical(TF.medication) == 'placebo');
tf = TF.occipital(rows);

% calculate fband limits
for a = 1:numel(TOI) - 1
    eval(TOI(a).definition);
    fband_lim(a) = window(2);  
end

% set limits of the figure
fig = figure(figure_counter);
loglog(x_log, closed_log, 'b:', 'LineWidth', 0.5)
yl = get(gca, 'ylim');
ylim([yl(1), yl(2) + 0.1])
xlim(window_visual)
hold on

% plot lower alpha window
r1 = rectangle('Position', [fband_lim(2), yl(1) + 0.01, fband_lim(4) - fband_lim(2), yl(2) + 0.1], 'FaceColor', [1 0.7 0.7], 'EdgeColor', 'none');
r2 = rectangle('Position', [fband_lim(4), yl(1) + 0.01, fband_lim(5) - fband_lim(4), yl(2) + 0.1], 'FaceColor', [1 0.6 0.6], 'EdgeColor', 'none');

% plot lines between frequency bands (= TOIs)
for a = 1:length(fband_lim)
    l(a) = line([fband_lim(a), fband_lim(a)], [yl(1), yl(2) + 0.1], 'Color', [0.6 0.6 0.6]); 
end
h = line([window_visual(1), window_visual(2)], [yl(1) + 0.01, yl(1) + 0.01], 'Color', [0.6 0.6 0.6]);

% plot data
pl(1) = loglog(x_log, open_log, 'Color', [0, 0, 0], 'LineWidth', 2.5, 'LineStyle', ':');
hold on
pl(2) = plot(x_log, closed_log, 'Color', [0, 0, 0], 'LineWidth', 2.5);
hold on

% plot IAF and TF
v(1) = line([iaf, iaf], [yl(1), yl(2) + 0.1], 'Color', colours2(4, :), 'LineWidth', 3);
v(2) = line([tf, tf], [yl(1), yl(2) + 0.1], 'Color', colours2(4, :), 'LineWidth', 3, 'LineStyle', '--');

% add names of frequency bands
% for a = 1:numel(TOI)
%     if a == 1
%         text(window_visual(1) + ((fband_lim(a) - 0.1)/2), yl(1) + 0.005, TOI(a).sign, 'fontsize', 14, 'color', [0, 0, 0]);
%     elseif a == 8
%         text((fband_lim(a-1) + ((window_visual(2) - fband_lim(a-1))/2)), yl(1) + 0.005, TOI(a).sign, 'fontsize', 14, 'color', [0, 0, 0]);
%     else
%         text((fband_lim(a-1) + ((fband_lim(a) - fband_lim(a-1))/2)), yl(1) + 0.005, TOI(a).sign, 'fontsize', 14, 'color', [0, 0, 0]);
%         hold on
%     end
% end

% add parameters  
xlabel('frequency (Hz)')
ylabel('relative amplitude')
set(gca, 'FontSize', 16)

% add legend 
legend(pl, {'eyes open', 'eyes closed'},...
    'FontSize', 16,'position', [0.65, 0.65, 0.22, 0.22], 'edgecolor', [0.6 0.6 0.6])
pause(5)

% save figure
figure_name = 'rsEEG_fbands_example';
savefig(figure_name)
saveas(fig, [figure_name '.png'])

clear data_visual_open open_log data_visual_closed closed_log x x_log yl fband_lim r1 r2 l h a pl v tf iaf figure_name

%% 6) extract mean amplitude values over frequency bands 
% prepare a long format table for R
fband_amplitude = table;
fband_amplitude.subject = (repelem(participant, length(medication) * length(time) * length(condition) * numel(TOI)))';
fband_amplitude.medication = repmat((repelem(medication, length(time) * length(condition) * numel(TOI)))', length(participant), 1);
fband_amplitude.time = repmat((repelem(time, length(condition) * numel(TOI)))', length(participant) * length(medication), 1);
fband_amplitude.condition = repmat((repelem(condition, numel(TOI)))', length(participant) * length(medication) * length(time), 1);
fband_amplitude.fband = repmat({TOI(:).band}', length(participant) * length(medication) * length(time) * length(condition), 1);

% calculate amplitudes of all subjects individually
% p = 1; m = 1; t = 1; c = 1; f = 1; 
for r = 1:numel(ROI)
    amplitude = [];
    for p = 1:length(participant)
        for m = 1:length(medication)
            for t = 1:length(time)
                for c = 1:length(condition)
                    for f = 1:numel(TOI)
                        % call individual fband data
                        window = cell2mat(fband{p, 2 + f}); 

                        % calculate amplitude
                        if f == 1 
                            % identify fband limits
                            x_start = ceil((window(1) - header_low.xstart)/header_low.xstep + 1);
                            x_end = ceil((window(2) - header_low.xstart)/header_low.xstep);

                            % perform the calculation on windowed data 
                            mean_amp = mean(squeeze(data_low(m, t, c, p, r, x_start : x_end)));

                            % append the result
                            amplitude = [amplitude; mean_amp];
                        else
                            % identify fband limits
                            x_start = ceil((window(1) - header_high.xstart)/header_high.xstep + 1);
                            x_end = ceil((window(2) - header_high.xstart)/header_high.xstep);

                            % perform the calculation on windowed data 
                            mean_amp = mean(squeeze(data_high(m, t, c, p, r, x_start : x_end)));
                            
                            % append the result
                            amplitude = [amplitude; mean_amp];
                        end
                    end
                end
            end
        end
    end
    statement = ['fband_amplitude.' ROI(r).area ' = amplitude;'];
    eval(statement)
end

% save the table
save('rsEEG_fband_amplitude.mat', 'fband_amplitude')
writetable(fband_amplitude, 'rsEEG_fband_amplitude.csv')

clear amplitude window x_start x_end mean_amp statement

%% 7) plot mean amplitudes 
% calculate group mean values
for m = 1:length(medication)
    for t = 1:length(time)
        for c = 1:length(condition)
            for r = 1:numel(ROI)
                for f = 1:numel(TOI)
                    avg_amp(m, t, c, r, f) = mean(fband_amplitude.(5 + r)(strcmp(fband_amplitude.medication, medication{m}) & strcmp(fband_amplitude.time, time{t}) &...
                        strcmp(fband_amplitude.condition, condition{c}) & strcmp(fband_amplitude.fband, TOI(f).band)));
                    avg_amp_std(m, t, c, r, f) = std(fband_amplitude.(5 + r)(strcmp(fband_amplitude.medication, medication{m}) & strcmp(fband_amplitude.time, time{t}) &...
                        strcmp(fband_amplitude.condition, condition{c}) & strcmp(fband_amplitude.fband, TOI(f).band)));
                    avg_amp_sem(m, t, c, r, f) = avg_amp_std(m, t, c, r, f)/sqrt(length(participant));
                end
            end
        end
    end
end
clear m t c r f 
disp(['Datasize: ' num2str(size(avg_amp))])

% barplots - plot by region
for r = 1:numel(ROI)
    for c = 1:length(condition)
        % choose data
        data_visual = [squeeze(avg_amp(1, 1, c, r, :))'; squeeze(avg_amp(1, 2, c, r, :))'; squeeze(avg_amp(2, 1, c, r, :))'; squeeze(avg_amp(2, 2, c, r, :))'];
        sem_visual = [squeeze(avg_amp_sem(1, 1, c, r, :))'; squeeze(avg_amp_sem(1, 2, c, r, :))'; squeeze(avg_amp_sem(2, 1, c, r, :))'; squeeze(avg_amp_sem(2, 2, c, r, :))'];
        
        % launch the figure
        fig = figure(figure_counter);
        hold on
        barplot = bar(data_visual', 'EdgeColor', 'none');
        for a = 1:size(data_visual, 1)
            barplot(a).FaceColor = colours2(a, :)
        end
        
        % plot errorbars
        ngroups = numel(TOI);
        nbars = size(data_visual, 1);
        groupwidth = min(0.8, nbars/(nbars + 1.5));
        for i = 1:nbars
            x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
            errorbar(x, data_visual(i, :), [], sem_visual(i, :), 'k', 'linestyle', 'none');
        end        
        
        % set other parameters
        fig_name = [ROI(r).area ' region, eyes ' condition{c}];
        title(fig_name, 'FontWeight', 'bold', 'FontSize', 18)
        legend(barplot, {'placebo - baseline', 'placebo - post med', 'alprazolam - baseline', 'alprazolam - post med'}, 'Location', 'northeast', 'fontsize', 16)
        xlabel('frequency band')
        ylabel('relative amplitude (\muV)')
        set(gca, 'xtick', 1:numel(TOI), 'xticklabel', {TOI.band})
        set(gca, 'Fontsize', 16)
        hold off
        
        % wait for size adjustment for PNG
        disp('Time to change the figure size, if necessary.')
        pause(5)

        % save figure
        fig_name = ['rsEEG_barplot_' ROI(r).area '_' condition{c}];
        savefig([fig_name '.fig'])
        saveas(fig, [fig_name '.png'])
        
        % update counter
        figure_counter = figure_counter + 1;
    end
end
clear r c fig fig_name ngroups nbars groupwidth fig barplot data_visual sem_visual

% normalize amplitude to baseline
for m = 1:length(medication)
    for c = 1:length(condition)
        for r = 1:numel(ROI)
            for f = 1:numel(TOI)
                for p = 1:length(participant)
                    amp_baseline = fband_amplitude.(5 + r)(strcmp(fband_amplitude.medication, medication{m}) & ... 
                        strcmp(fband_amplitude.time, time{1}) & ... % choose time pre medication
                        strcmp(fband_amplitude.condition, condition{c}) & ...
                        strcmp(fband_amplitude.fband, TOI(f).band) & ...
                        fband_amplitude.subject == participant(p));
                    amp_raw = fband_amplitude.(5 + r)(strcmp(fband_amplitude.medication, medication{m}) & ... 
                        strcmp(fband_amplitude.time, time{2}) & ... % choose time post medication
                        strcmp(fband_amplitude.condition, condition{c}) & ...
                        strcmp(fband_amplitude.fband, TOI(f).band) & ...
                        fband_amplitude.subject == participant(p));
                    amp_norm(m, c, r, f, p) = amp_raw / amp_baseline * 100;
                end
            end
        end
    end
end
clear m c r f p amp_baseline amp_raw

% plot relative amplitude by fband
for f = 1:numel(TOI)
    for c = 1:length(condition)
        % choose data
        data_visual = squeeze(amp_norm(:, c, :, f, :));
       
        % launch the figure
        fig = figure(figure_counter);
        hold on
        
        % add no-change line
        xlim([0.75 size(data_visual, 2) + 0.25])
        line([0.75 size(data_visual, 2) + 0.25], [100, 100], 'Color', [0.75 0.75 0.75], ...
            'LineStyle', '--', 'LineWidth', 2);
                              
        % plot data with errorbars
        x = 1:size(data_visual, 2);    
        for m = 1:size(data_visual, 1)  % medication
            % calculate the data and 95% CI
            for r = 1:size(data_visual, 2); % ROIs  
                y(r) = mean(squeeze(data_visual(m, r, :)));
                CI(r) = std(squeeze(data_visual(m, r, :))) / sqrt(length(participant)) * z;
            end

            % plot
            err(m) = errorbar(x, y, CI)

            % add parameters
            err(m).Color = colours2((m - 1)*2 + 2, :);
            err(m).LineWidth = 1.2;
            err(m).Marker = 'o';
            err(m).MarkerFaceColor = colours2((m - 1)*2 + 2, :);
            err(m).MarkerSize = 10;
        end
        
        % set other parameters
        fig_name = [TOI(f).band ' frequency band, eyes ' condition{c}];
        title(fig_name, 'FontWeight', 'bold', 'FontSize', 16)
        set(gca, 'Fontsize', 14)
        legend(err, medication, 'Location', 'northeast', 'fontsize', 14, 'EdgeColor', 'none')
        xlabel('region of interest')
        set(gca, 'xtick', 1:size(data_visual, 2), 'xticklabel', {'frontal' 'central' 'left' 'right' 'occipital'})
        ylabel('relative amplitude (% baseline)')
        hold off
        
        % wait for size adjustment for PNG
        disp('Time to change the figure size, if necessary.')
        pause(2)

        % save figure
        fig_name = ['rsEEG_change_' TOI(f).band '_' condition{c}];
        savefig([fig_name '.fig'])
        saveas(fig, [fig_name '.png'])
        
        % update counter
        figure_counter = figure_counter + 1;
    end
end
clear f c m r x y CI err fig fig_name data_visual sem_visual amp_norm

%% 8) beta band 
% choose data
beta_fbands = {'beta1' 'beta2'};
for m = 1:length(medication)
    for t = 1:length(time)
        for c = 1:length(condition)
            for r = 1:numel(ROI)
                for b = 1:length(beta_fbands)
                    % extract amplitude of individual beta sub-bands
                    rows = (categorical(fband_amplitude.medication) == medication{m} & ...
                        categorical(fband_amplitude.time) == time{t} & ...
                        categorical(fband_amplitude.condition) == condition{c} & ...
                        categorical(fband_amplitude.fband) == beta_fbands{b});
                    data_beta(m, t, c, r, b, :) = fband_amplitude{rows, ROI(r).area};
                end
            
                % calculate broad band beta amplitude
                for p = 1:length(participant)
                    data_beta(m, t, c, r, b+1, p) =  mean(data_beta(m, t, c, r, [1:b], p));
                end
            end
        end
    end 
end
clear m t c r b rows
disp(['Datasize: ' num2str(size(data_beta))])

% extract amplitude from all beta bands --> central region, eyes open
beta_fbands{end + 1} = 'broad';
for m = 1:length(medication)
    % raw beta
    for t = 1:length(time)
        for p = 1:length(participant)
            for b = 1:length(beta_fbands)
                beta(m, t, p, b) = squeeze(data_beta(m, t, 1, 2, b, p));
            end
        end
    end
    
    % change --> normalized beta
    for p = 1:length(participant)
        for b = 1:length(beta_fbands)
            beta_change(m, p, b) = (beta(m, 2, p, b)/beta(m, 1, p, b))*100 - 100;
        end
    end
end
save('beta.mat', 'beta'); save('beta_change.mat', 'beta_change'); 
clear m t p b

% plot relative amplitude by fband
col = colours([2 4], :);
beta_fbands_range = {'low beta = sigma peak' 'high beta' 'broad band beta'};
for b = 1:length(beta_fbands)
    % choose the data
    data_visual = [];
    for m = 1:length(medication)
        data_i = squeeze(beta_change(m, :, b))';
        data_visual = cat(2, data_visual, data_i);
    end
    
    % launch the figure
    fig = figure(figure_counter);
    hold on

    % plot group boxplot
    boxplot(data_visual, 'color', col)
    
    % add zero line
    xl = get(gca, 'xlim');
    h = line([xl(1) xl(2)], [0, 0], 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 2);

    % plot the lines
    for p = 1:length(participant)
        p_line(p) = plot([1 2], data_visual(p, [1 2]), '-o',...
            'Color', [0.75, 0.75, 0.75],...
            'MarkerSize', 10,...
            'MArkerEdge', 'none');
        hold on
    end

    % plot the markers
    for s = 1:size(data_visual, 2)
        scat(s) = scatter(repelem(s, length(participant)), data_visual(:, s),...
            75, col(s, :), 'filled');
        hold on
    end

    % add parameters
    figure_title = ['BETA POWER CHANGE: ' beta_fbands_range{b}];
    figure_name = ['rsEEG_beta_change_' beta_fbands{b}];
    yl = get(gca, 'ylim');
    set(gca, 'xtick', [1 2], 'xticklabel', {'placebo' 'alprazolam'})
    set(gca, 'Fontsize', 16)
    title(figure_title, 'FontWeight', 'bold', 'FontSize', 18)
    xlabel('medication'); ylabel('power change (% baseline)');
    ylim([yl(1), yl(2)])
    hold off

    % save figure
    savefig([figure_name '.fig'])
    saveas(fig, [figure_name '.png'])

    % update counter
    figure_counter = figure_counter + 1;
    clear b m fig data_i data_visual figure_title figure_name yl xl p p_line h s scat 
end
clear col beta_fbands_range 

%% 8) delta band 
% choose data
for m = 1:length(medication)
    for t = 1:length(time)
        for c = 1:length(condition)
            for r = 1:numel(ROI)
                % extract amplitude of delta band
                rows = (categorical(fband_amplitude.medication) == medication{m} & ...
                    categorical(fband_amplitude.time) == time{t} & ...
                    categorical(fband_amplitude.condition) == condition{c} & ...
                    categorical(fband_amplitude.fband) == 'delta');
                data_delta(m, t, c, r, :) = fband_amplitude{rows, ROI(r).area};
            end
        end
    end 
end
clear m t c r rows
disp(['Datasize: ' num2str(size(data_delta))])

% extract delta amplitude --> frontal region, eyes open
for m = 1:length(medication)
    % raw delta
    for t = 1:length(time)
        for p = 1:length(participant)
            delta(m, t, p) = squeeze(data_delta(m, t, 1, 1, p));
        end
    end
    
    % change --> normalized beta
    for p = 1:length(participant)
        delta_change(m, p) = (delta(m, 2, p)/delta(m, 1, p))*100 - 100;
    end
end
save('delta.mat', 'delta'); save('delta_change.mat', 'delta_change'); 
clear m t p 

% plot relative amplitude 
% choose the data
data_visual = [];
for m = 1:length(medication)
    data_i = squeeze(delta_change(m, :))';
    data_visual = cat(2, data_visual, data_i);
end

% launch the figure
fig = figure(figure_counter);
hold on

% plot group boxplot
col = colours([2 4], :);
boxplot(data_visual, 'color', col)

% add zero line
xl = get(gca, 'xlim');
h = line([xl(1) xl(2)], [0, 0], 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 2);

% plot the lines
for p = 1:length(participant)
    p_line(p) = plot([1 2], data_visual(p, [1 2]), '-o',...
        'Color', [0.75, 0.75, 0.75],...
        'MarkerSize', 10,...
        'MArkerEdge', 'none');
    hold on
end

% plot the markers
for s = 1:size(data_visual, 2)
    scat(s) = scatter(repelem(s, length(participant)), data_visual(:, s),...
        75, col(s, :), 'filled');
    hold on
end

% add parameters
figure_title = 'DELTA POWER CHANGE';
figure_name = 'rsEEG_delta_change';
yl = get(gca, 'ylim');
set(gca, 'xtick', [1 2], 'xticklabel', {'placebo' 'alprazolam'})
set(gca, 'Fontsize', 16)
title(figure_title, 'FontWeight', 'bold', 'FontSize', 18)
xlabel('medication'); ylabel('power change (% baseline)');
ylim([yl(1), yl(2)])
hold off

% save figure
savefig([figure_name '.fig'])
saveas(fig, [figure_name '.png'])

% update counter
figure_counter = figure_counter + 1;
clear col fig data_i data_visual figure_title figure_name yl xl p p_line h s scat  

%% 10) alpha attenuation coeficient  
% ----- extract AAC -----
% choose data - individual alpha subbands + broad alpha band
alpha_fbands = {'alpha1' 'alpha2' 'alpha3'};
for m = 1:length(medication)
    for t = 1:length(time)
        for c = 1:length(condition)
            for a = 1:length(alpha_fbands)
                % extract amplitude of individual alpha sub-bands
                rows = (categorical(fband_amplitude.medication) == medication{m} & ...
                    categorical(fband_amplitude.time) == time{t} & ...
                    categorical(fband_amplitude.condition) == condition{c} & ...
                    categorical(fband_amplitude.fband) == alpha_fbands{a});
                data_aac(m, t, c, a, :) = fband_amplitude{rows, 'occipital'};
            end
            
            % calculate broad band alpha amplitude
            for p = 1:length(participant)
                data_aac(m, t, c, a+1, p) =  mean(data_aac(m, t, c, [1:a], p));
            end
        end
    end 
end
% update alpha bands
alpha_fbands{end + 1} = 'broad';
alpha_fbands_txt = {'low alpha 1' 'low alpha 2'  'high alpha 3' 'broad alpha'};

% calculate ACC
for m = 1:length(medication)
    % calculate individual AAC
    for t = 1:length(time)
        for a = 1:length(alpha_fbands)
            for p = 1:length(participant)
                AAC(m, t, a, p) = squeeze(data_aac(m, t, 2, a, p)) / squeeze(data_aac(m, t, 1, a, p));
            end
        end
    end
    
    % calculate AAC change 
    for a = 1:length(alpha_fbands)
        for p = 1:length(participant)
            AAC_change(m, a, p) = (AAC(m, 2, a, p)/AAC(m, 1, a, p))*100 - 100;
        end
    end
end
save('AAC.mat', 'AAC'); save('AAC_change.mat', 'AAC_change'); 

% ----- plot box + scatter plot -----
% individual ACC - for each alpha band separately
for a = 1:length(alpha_fbands)
    % choose the data
    data_visual = [];
    for m = 1:length(medication)
        for t = 1:length(time)
            data_i = squeeze(AAC(m, t, a, :));
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
    figure_title = ['ALPHA ATTENUATION COEFFICIENT: ' alpha_fbands_txt{a}];
    figure_name = ['rsEEG_ACC_' alpha_fbands{a}];
    yl = get(gca, 'ylim');
    set(gca, 'xtick', 1:4, 'xticklabel', {'pre' 'post' 'pre' 'post'})
    set(gca, 'Fontsize', 16)
    title(figure_title, 'FontWeight', 'bold', 'FontSize', 18)
    xlabel('time relative to medication'); ylabel('alpha attenuation coefficient');
    ylim([0, yl(2) + 1])
    hold on

    % add text
    txt(1) = text(1.3, yl(2) + 0.5, 'placebo', 'fontsize', 16, 'color', colours2(2, :));
    txt(2) = text(3.3, yl(2) + 0.5, 'alprazolam', 'fontsize', 16, 'color', colours2(4, :));
    hold off

    % save the figure
    pause(5)        
    savefig([figure_name '.fig'])
    saveas(fig, [figure_name '.png'])

    % update the counter
    figure_counter = figure_counter + 1;        
end

% ACC change - for each alpha band separately
for a = 1:length(alpha_fbands)
    % choose the data
    data_visual = [];
    for m = 1:length(medication)
        data_i = squeeze(AAC_change(m, a, :));
        data_visual = cat(2, data_visual, data_i);
    end

    % plot group boxplot
    col = colours([2 4], :);
    fig = figure(figure_counter);        
    boxplot(data_visual, 'color', col)
    hold on
    
    % add zero line
    xl = get(gca, 'xlim');
    h = line([xl(1) xl(2)], [0, 0], 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 2);
    hold on

    % plot the lines
    for p = 1:length(participant)
        p_open(p) = plot([1 2], data_visual(p, [1 2]), '-o',...
            'Color', [0.75, 0.75, 0.75],...
            'MarkerSize', 10,...
            'MArkerEdge', 'none');
        hold on
    end

    % plot the markers
    for b = 1:size(data_visual, 2)
        scat(b) = scatter(repelem(b, length(participant)), data_visual(:, b),...
            75, col(b, :), 'filled');
        hold on
    end

    % add parameters
    figure_title = ['AAC CHANGE: ' alpha_fbands_txt{a}];
    figure_name = ['rsEEG_ACCchange_' alpha_fbands{a}];
    yl = get(gca, 'ylim');
    set(gca, 'xtick', [1 2], 'xticklabel', {'placebo' 'alprazolam'})
    set(gca, 'Fontsize', 16)
    title(figure_title, 'FontWeight', 'bold', 'FontSize', 18)
    xlabel('medication'); ylabel('AAC change');
    ylim([yl(1), yl(2)])
    hold off

    % save the figure
    pause(5)        
    savefig([figure_name '.fig'])
    saveas(fig, [figure_name '.png'])

    % update the counter
    figure_counter = figure_counter + 1;        
end

%% SPECTRAL EXPONENT: parameters
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
