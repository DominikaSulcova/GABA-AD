%% RS-EEG - EXTRACTION OF MEAN FREQUENCY BAND AMPLITUDE 
% Written by Dominika for GABA-AD project (2021)
% 
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
% 7) Plots barplots with mean amplitudes for each ROI
% 8) Calculates alpha attenuation coefficient (AAC)
%       - for individual alpha subbands + broad alpha band
%       - calculates individual AAC and AAC change 
%       - plots box + scatter plots, saves figures
%       --> saves values as 'AAC.mat' and 'AAC.csv'

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
load('colours.mat'); load('colours2.mat'); 

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
            ylabel('relative amplitude (µV)')
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
disp(['Datasize: ' num2str(size(avg_amp))])

% plot by region
for r = 1:numel(ROI)
    for c = 1:length(condition)
        % choose data
        data_visual = [squeeze(avg_amp(1, 1, c, r, :))'; squeeze(avg_amp(1, 2, c, r, :))'; squeeze(avg_amp(2, 1, c, r, :))'; squeeze(avg_amp(2, 2, c, r, :))'];
        sem_visual = [squeeze(avg_amp_sem(1, 1, c, r, :))'; squeeze(avg_amp_sem(1, 2, c, r, :))'; squeeze(avg_amp_sem(2, 1, c, r, :))'; squeeze(avg_amp_sem(2, 2, c, r, :))'];
        
        % launch the figure
        fig = figure(figure_counter);
        barplot = bar(data_visual', 'EdgeColor', 'none');
        colormap(colours2)
        hold on
        
        % plot errorbars
        ngroups = numel(TOI);
        nbars = size(data_visual, 1);
        groupwidth = min(0.8, nbars/(nbars + 1.5));
        for i = 1:nbars
            x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
            errorbar(x, data_visual(i, :), [], sem_visual(i, :), 'k', 'linestyle', 'none');
            hold on
        end        
        
        % set other parameters
        fig_name = [ROI(r).area ' region, eyes ' condition{c}];
        title(fig_name, 'FontWeight', 'bold', 'FontSize', 18)
        legend(barplot, {'placebo - baseline', 'placebo - post med', 'alprazolam - baseline', 'alprazolam - post med'}, 'Location', 'northeast', 'fontsize', 16)
        xlabel('frequency band')
        ylabel('relative amplitude (µV)')
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

clear data_visual sem_visual

%% 8) alpha attenuation coeficient  
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
                data_aac(m, t, c, a+1, p) =  sum(data_aac(m, t, c, [1:a], p));
            end
        end
    end 
end
% update alpha bands
alpha_fbands = [alpha_fbands, {'broad'}];
alpha_fbands_txt = {'low alpha 1' 'low alpha 2'  'high alpha 3' 'broad alpha'};

% calculate ACC
for m = 1:length(medication)
    % calculate individual AAC
    for t = 1:length(time)
        for a = 1:length(alpha_fbands)
            for p = 1:length(participant)
                ACC(m, t, a, p) = squeeze(data_aac(m, t, 2, a, p)) / squeeze(data_aac(m, t, 1, a, p));
            end
        end
    end
    
    % calculate AAC change 
    for a = 1:length(alpha_fbands)
        for p = 1:length(participant)
            ACC_change(m, a, p) = (ACC(m, 2, a, p)/ACC(m, 1, a, p))*100 - 100;
        end
    end
end
save('ACC.mat', 'ACC'); save('ACC_change.mat', 'ACC_change'); 

% ----- plot box + scatter plot -----
% individual ACC - for each alpha band separately
for a = 1:length(alpha_fbands)
    % choose the data
    data_visual = [];
    for m = 1:length(medication)
        for t = 1:length(time)
            data_i = squeeze(ACC(m, t, a, :));
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
        data_i = squeeze(ACC_change(m, a, :));
        data_visual = cat(2, data_visual, data_i);
    end

    % plot group boxplot
    col = colours2([2 4], :);
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

%% 9) beta increase  
% ----- extract BI -----
% choose data - individual beta subbands + broad beta band
beta_fbands = {'beta1' 'beta2'};
for m = 1:length(medication)
    for t = 1:length(time)
        for c = 1:length(condition)
            for b = 1:length(beta_fbands)
                % extract amplitude of individual alpha sub-bands
                rows = (categorical(fband_amplitude.medication) == medication{m} & ...
                    categorical(fband_amplitude.time) == time{t} & ...
                    categorical(fband_amplitude.condition) == condition{c} & ...
                    categorical(fband_amplitude.fband) == beta_fbands{b});
                data_bi(m, t, c, b, :) = fband_amplitude{rows, 'occipital'};
            end
            
            % calculate broad band alpha amplitude
            for p = 1:length(participant)
                data_bi(m, t, c, b+1, p) =  sum(data_bi(m, t, c, [1:b], p));
            end
        end
    end 
end
% update alpha bands
beta_fbands = [beta_fbands, {'broad'}];
beta_fbands_txt = {'low beta 1' 'high beta 2' 'broad beta'};

% calculate BI
for m = 1:length(medication)
    % calculate individual BI
    for c = 1:length(condition)
        for b = 1:length(beta_fbands)
            for p = 1:length(participant)
                BI(m, c, b, p) = squeeze(data_bi(m, 2, c, b, p)) / squeeze(data_bi(m, 1, c, b, p)) * 100;
            end
        end
    end
end
save('BI.mat', 'BI')

% ----- plot box + scatter plot -----
% BI - for each beta band separately
for a = 1:length(beta_fbands)
    for c = 1:length(condition)
        % choose the data
        data_visual = [];
        for m = 1:length(medication)
            data_i = squeeze(BI(m, c, a, :));
            data_visual = cat(2, data_visual, data_i);
        end

        % plot group boxplot
        col = colours2([2 4], :);
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
        figure_title = ['Beta increase - eyes ' condition{c} ' - ' beta_fbands_txt{a}];
        figure_name = ['rsEEG_BI_' condition{c} '_' beta_fbands{a}];
        yl = get(gca, 'ylim');
        set(gca, 'xtick', [1 2], 'xticklabel', {'placebo' 'alprazolam'})
        set(gca, 'Fontsize', 16)
        title(figure_title, 'FontWeight', 'bold', 'FontSize', 18)
        xlabel('medication'); ylabel('amplitude increase (% baseline)');
        ylim([yl(1), yl(2)])
        hold off

        % save the figure
        pause(5)        
        savefig([figure_name '.fig'])
        saveas(fig, [figure_name '.png'])

        % update the counter
        figure_counter = figure_counter + 1;      
    end
end
