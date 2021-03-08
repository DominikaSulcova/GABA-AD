%% MEP PROCESSING
% Written by Dominika for the GABA-AD project (2020 - 2021)
% 
% ----- Discards MEP epochs with baseline activity ----- 
% 1) Loops through the datset, in one loop:
%       - calculates average RMS of the baseline  
%       - discards all epochs that have baseline RMS larger than 
%         average RMS + <allow_sd> * sd
%    In case that there are discarded epochs, recalculates --> next cycle
% 2) Discards all epochs that at any point depass threshold value 
% 3) Saves the dataset, updates outcome table + saves outcome figure
% 
% ----- Extracts peak-to-peak amplitude ----- 
% 4) Calculates the amplitude for each epoch
% 5) Identifies zero epochs - p2p amplitude < 2 * threshold
%       - splits in two variables with and without zero epochs  
%       - mean p2p amplitude --> in the 'output' table
%       - saves sorted data as LW dataset --> prefix 'nozero'
% 
% ----- Group analysis ----- 
% 6) Unblinds medication for sessions
% 7) Plots group values 
%       - drug-induced change in spTMS MEPs --> box + scatterplots 
%       - SICI

%% parameters
clear all
clc

% dataset 
participant = [1:20];
session = {'S1' 'S2'};
medication = {'placebo' 'alprazolam'};
time = {'pre' 'post'};
stimulus = {'TS' 'ppTMS'};

% output
prefix_1 = 'clean';
prefix_2 = 'nozero';
output_name = 'GABA_MEP';

% baseline filter
baseline = -0.1;
threshold = 15;
allow_sd = 3;

% amplitude
window = [0.015 0.05];

% visualization
figure_counter = 1;
load('colours2.mat');

%% BASELINE ACTIVITY
% prepare an output table
if exist([output_name '.mat']) == 0
    output = table; 
    output.subject = zeros(0); output.session = zeros(0); output.time = zeros(0); output.stimulus = zeros(0);
    output.discarded = zeros(0); output.percent = zeros(0); output.cycles = zeros(0);
    output.threshold = zeros(0);
else
    load([output_name '.mat'])
end
    
% loop through participants, sessions and timepoints
for p = 1:length(participant)
    for n = 1:length(session)
        for t = 1:length(time)
            for s = 1:length(stimulus)
                % load header and dataset 
                dataset_name = ['complete ' stimulus{s} ' YC' num2str(participant(p)) ' ' session{n} ' ' time{t} ' EMG'];
                load([dataset_name '.lw6'],'-mat');      
                load([dataset_name '.mat']); 

                % prepare baseline data for visualization, detrend
                data_visual = squeeze(data(:, 1, 1, 1, 1, 1 : ceil(-baseline / header.xstep)));
                for a = 1:size(data_visual, 1)
                    data_visual(a, :) = detrend(data_visual(a, :));
                end
                data_visual_orig = data_visual;

                % prepare x 
                x = baseline : header.xstep : 0; 

                % outcome vector
                discarded = [];

                %% 1) run automatic RMS + SD removal 
                go = 1; cycle = 1; 
                while go 
                    discarded_pos = [];

                    % calculate average RMS values
                    avg_rms = mean(rms(data_visual'));
                    avg_sd = std(rms(data_visual'));
                    cutoff = avg_rms + allow_sd * avg_sd;

                    % plot original dataset
                    fig = figure(figure_counter); 
                    set(gcf,'units','normalized','outerposition',[0 0 1 1])
                    subplot(3, 2, [1 2])
                    plot(x, data_visual, 'color', [0.45, 0.45, 0.45], 'linewidth', 1.5)                
                    suptitle(['Subject ' num2str(participant(p)) ' - ' session{n} ' - ' time{t} ' - ' stimulus{s} ' : CYCLE ' num2str(cycle)])
                    title(['Original dataset - average RMS ' num2str(avg_rms) ' \muV'])
                    set(gca,'Fontsize',16); ylabel('amplitude (\muV)');
                    hold on

                    % loop through trials
                    for e = 1:length(header.events)
                        % check if the rms of current event fits into the limis 
                            current_rms = rms(data_visual(e, :)');
                            if current_rms > cutoff
                                discarded_pos = [discarded_pos e];
                                discarded = [discarded header.events(e).epoch];
                                ditch = true;
                            else
                                ditch = false;
                            end                        

                        % plot the event to the appropriate axes                
                        if ditch   
                            subplot(3, 2, [5 6])
                            plot(x, data_visual(e, :), 'Color', [0.99, 0.3, 0.2], 'LineWidth', 1.5)
                            hold on
                        else
                            subplot(3, 2, [3 4])
                            plot(x, data_visual(e, :), 'Color', [0, 0, 0], 'LineWidth', 1.5)
                            hold on
                        end
                    end

                    % add parameters to axes
                    subplot(3, 2, [3 4])
                    title(['Kept epochs: ' num2str(length(header.events) - length(discarded_pos))])
                    set(gca,'Fontsize',16); ylabel('amplitude (\muV)');
                    hold on

                    subplot(3, 2, [5 6])
                    title(['Discarded epochs: ' num2str(length(discarded_pos))])
                    set(gca,'Fontsize',16); xlabel('time (s)'); ylabel('amplitude (\muV)');
                    hold on

                    % ------------------ cycle automatically to 0 discarded epochs ------------------
                    if ~isempty(discarded_pos)
                        % remove indicated epochs from data, update header
                        data(discarded_pos, :, :, :, :, :) = [];
                        header.datasize(1) = header.datasize(1) - length(discarded_pos);
                        header.events(discarded_pos)= [];

                        % remove indicated epochs from visual dataset for
                        % future filtration cycles
                        data_visual(discarded_pos , :) = [];

                        % continue the cycle
                        pause(1); clf;
                        cycle = cycle + 1;  
                    else
                        % close current figure
                        pause(1); close(fig);   

                        % plot original dataset
                        fig = figure(figure_counter); 
                        set(gcf,'units','normalized','outerposition',[0 0 1 1])
                        subplot(3, 2, 3)
                        plot(x, data_visual_orig, 'color', [0, 0, 0], 'linewidth', 1.5)                       
                        title(['Original dataset: ' num2str(size(data_visual_orig, 1)) ' epochs'])
                        set(gca,'Fontsize',14); ylabel('amplitude (\muV)');
                        hold on  

                        % plot filtered dataset
                        subplot(3, 2, 5)
                        plot(x, data_visual, 'Color', [0, 0, 0], 'linewidth', 1.5)
                        xlim = get(gca,'xlim');                        
                        title(['RMS + SD: ' num2str(length(header.events)) ' epochs kept, ' num2str(cycle - 1) ' cycles performed'])
                        set(gca,'Fontsize',14); ylabel('amplitude (\muV)'); xlabel('time (s)');
                        hold on

                        % add threshold 
                        subplot(3, 2, 5)
                        l(1) = line([baseline, 0], [threshold, threshold], 'LineWidth', 1.5, 'Color', [0.99, 0.3, 0.2], 'LineStyle', '--');
                        l(2) = line([baseline, 0], [-threshold, -threshold], 'LineWidth', 1.5, 'Color', [0.99, 0.3, 0.2], 'LineStyle', '--');
                        text(xlim(1) + 0.005 , - threshold + 4 ,['threshold = ' num2str(threshold) ' µV'], 'Fontsize', 14, 'color', [0.99, 0.3, 0.2])
                        hold on

                        % plot discarded epochs - RMS + SD
                        if ~isempty(discarded)
                            subplot(3, 2, 4)
                            plot(x, data_visual_orig(discarded, :), 'Color', [0.99, 0.3, 0.2], 'LineWidth', 1.5) 
                            title(['RMS + SD: ' num2str(length(discarded)) ' epochs discarded'])
                            set(gca,'Fontsize',14)
                            hold on
                        else
                            subplot(3, 2, 4)
                            title('RMS + SD: No epoch discarded'); set(gca,'Fontsize',14)
                            hold on
                        end                     

                        % exit the while loop
                        go = 0;
                    end

                    % --------------------- choose final cycles manually ---------------------
    %                     answer = questdlg('Do you want to repeat the cycle?', 'Cycles',...
    %                     'Yes, keep them comin', 'No, save and continue to next dataset', 'No, save and continue to next dataset');
    %                     switch answer
    %                         case 'Yes, keep them comin'
    %                             cycle = cycle + 1; 
    %                             clf
    %                         case 'No, save and continue to next dataset' 
    %                             % close currentfigure
    %                             close(fig)     
    %                             
    %                          % plot the final figure and wait 5s
    %                         fig = figure(figure_counter); 
    %                         set(gcf,'units','normalized','outerposition',[0 0 1 1])
    %                         subplot(3, 2, [1 2])
    %                         plot(x, data_visual_orig)
    %                         h = suptitle(['Subject ' participant{p} ' - ' session{s} ' - ' time{t} ' : ' num2str(cycle - 1) ' cycles performed'])
    %                         set(h,'FontSize',20)
    %                         title(['Original dataset - ' num2str(size(data_visual_orig, 1)) ' epochs'])
    %                         set(gca,'Fontsize',16)
    %                         hold on  
    % 
    %                         subplot(3, 2, [3 4])
    %                         plot(x, data_visual, 'Color', [0, 0, 0], 'LineWidth', 1)
    %                         ylim = get(gca,'ylim'); xlim = get(gca,'xlim');                        
    %                         text(xlim(1) + 0.005 , ylim(2) - 2 ,['Final average RMS: ' num2str(avg_rms) ' µV'], 'Fontsize', 16)
    %                         title(['Kept epochs: ' num2str(length(header.events))])
    %                         set(gca,'Fontsize',16)
    %                         hold on
    % 
    %                         if ~isempty(discarded)
    %                             subplot(3, 2, [5 6])
    %                             plot(x, data_visual_orig(discarded, :), 'Color', [0.99, 0.3, 0.2], 'LineWidth', 1.5)
    %                             title(['Discarded epochs: ' num2str(length(discarded))])    
    %                             set(gca,'Fontsize',16)
    %                             hold off
    %                         end 
    %                             
    %                             pause(5)
    %                             
    %                             % save and close the figure
    %                             figure_name = ['filtered_' method{mode} '_' participant{p} '_' session{s} '_' time{t}]
    %                             savefig(fig, [figure_name '.fig'])
    %                             close(fig)
    %                             
    %                             % update the counter
    %                             figure_counter = figure_counter+ 1;
    %                             
    %                             % exit the while loop
    %                             go = 0;
    %                     end
    %                     % -----------------------------------------------------------------------
                end

                %% 2) remove epochs that depass the threshold
                % loop through left trials
                for e = 1:length(header.events)
                    % check if tthe maximum value across baseline datapoints
                    % fits under the threshold
                        current_max = max(abs(data_visual(e, :)));
                        if current_max > threshold
                            discarded_pos = [discarded_pos e];
                            discarded = [discarded header.events(e).epoch];
                            ditch = true;
                        else
                            ditch = false;
                        end                        

                    % plot the event to the appropriate axes                
                    if ~ditch   
                        subplot(3, 2, [1 2])
                        plot(x, data_visual(e, :), 'Color', [0, 0.45, 0.74], 'LineWidth', 1.5)
                        hold on
                    else
                        subplot(3, 2, 6)
                        plot(x, data_visual(e, :), 'Color', [0.99, 0.3, 0.2], 'LineWidth', 1.5)
                        hold on
                    end
                end

                % remove indicated epochs from visual dataset
                data_visual(discarded_pos , :) = [];

                % add parameters to axes
                final_rms = mean(rms(data_visual'));
                subplot(3, 2, [1 2])
                set(gca,'Fontsize',14); ylabel('amplitude (\muV)');  
                title(['Subject ' num2str(participant(p)) ', ' session{n} ', ' time{t} ', ' stimulus{s} ' - FINAL DATASET : ' num2str(length(header.events) - length(discarded_pos)) ' epochs kept, ' ...
                    num2str(length(discarded)) ' discarded - final average RMS ' num2str(final_rms) ' \muV'], 'Fontsize', 18)            
                hold on

                subplot(3, 2, 6)
                title(['THRESHOLD : ' num2str(length(discarded_pos)) '  epochs discarded'])
                set(gca,'Fontsize',14); xlabel('time (s)');
                hold on

                pause(3)

                %% 3) save outcome variables           
                % number final epochs
                for i = 1:header.datasize(1)
                    header.events(i).epoch = i;
                end

                % modify and save header
                header.datasize(1) = header.datasize(1) - length(discarded_pos);
                header.events(discarded_pos)= [];
                header.events = header.events';
                header.name = [prefix_1 ' ' dataset_name];
                save([header.name '.lw6'], 'header')         

                % modify and save data
                data(discarded_pos, :, :, :, :, :) = [];
                save([header.name '.mat'], 'data')

                % fill in the outcome table and save
                M = table;
                M.subject = participant(p); M.session = session(n); M.time = time(t); M.stimulus = stimulus(s);
                M.discarded = {sort(discarded)}; 
                M.percent = round(100 - (size(data_visual, 1) / size(data_visual_orig, 1)) * 100, 2); 
                M.cycles = cycle - 1; 
                M.threshold = threshold;

                output = [output; M];

                save([output_name '.mat'], 'output') 

                % save and close the figure
                figure_name = [prefix_1 '_YC' num2str(participant(p)) '_' session{n} '_' time{t} '_' stimulus{s}];
                savefig(fig, [figure_name '.fig'])
                close(fig)

                % update the counter
                figure_counter = figure_counter + 1;
            end
        end
    end
end
xstep = header.xstep;
clear avg_rms avg_sd current_max current_rms cutoff cycle data data_visual data_visual_orig dataset_name ...
    discarded discarded_pos ditch fig figure_name final_rms go header M xlim l x
clear a e i n p s t 

%% PEAK-TO-PEAK AMPLITUDE
% prepare a table for amplitude results
amplitudes = table; 
amplitudes.subject = zeros(0); amplitudes.session = zeros(0); amplitudes.time = zeros(0); amplitudes.stimulus = zeros(0);
amplitudes.zero = zeros(0); amplitudes.amp_zero = zeros(0); amplitudes.nozero = zeros(0); amplitudes.amp_nozero = zeros(0);

% loop through datasets 
for p = 1:length(participant)
    for n = 1:length(session)
        for t = 1:length(time)             
            for s = 1:length(stimulus)
                %% 4) extract p2p amplitude
                % load data and header
                load([prefix_1 ' complete ' stimulus{s} ' YC' num2str(participant(p)) ' ' session{n} ' ' time{t} ' EMG.mat']);
                load([prefix_1 ' complete ' stimulus{s} ' YC' num2str(participant(p)) ' ' session{n} ' ' time{t} ' EMG.lw6'], '-mat');
                
                % loop through epochs
                for e = 1:size(data, 1)
                    % choose the window
                    x_start = ceil((window(1) - header.xstart)/header.xstep);
                    x_end = ceil((window(2) - header.xstart)/header.xstep);
                    
                    % identify extremes
                    y_max(e) = max(squeeze(data(e, 1, 1, 1, 1, x_start:x_end)));
                    y_min(e) = min(squeeze(data(e, 1, 1, 1, 1, x_start:x_end)));
                    
                    % calculate amplitude 
                    amp(e) = y_max(e) - y_min(e); 
                end                
                
                %% 5) identify zero response epochs                 
                % loop through epochs
                for e = 1:length(amp)
                    if amp(e) > 2 * threshold
                        zero_index(e) = true; 
                    else
                        zero_index(e) = false; 
                    end
                end
                
                % create new dataset with no zero epochs
                data_nozero = data(zero_index, :, :, :, :, :); 
                
                % extract and append mean values
                A = table;
                A.subject = participant(p); A.session = session(n); A.time = time(t); A.stimulus = stimulus(s);
                A.zero = size(data, 1); A.amp_zero = mean(amp); 
                A.nozero = size(data_nozero, 1); A.amp_nozero = mean(amp(zero_index)); 
                
                amplitudes = [amplitudes; A];
                
                % adjust and save header
                header.name = [prefix_2 ' ' header.name];
                header.datasize(1) = size(data_nozero, 1);
                header.events = header.events(1:header.datasize(1));
                save([header.name '.lw6'], 'header')
                
                % save new dataset 
                data = data_nozero;
                save([header.name '.mat'], 'data')
                
                clear y_max y_min amp zero_index            
            end 
        end
    end
end
clear A data header data_nozero x_start x_end 
clear p n t s e

% make sure participants are ordered in result tables
output = sortrows(output, [1:4]);
amplitudes = sortrows(amplitudes, [1:4]);

% append results
output.epochs_zero = amplitudes.zero;
output.amplitude_zero = amplitudes.amp_zero;
output.epochs_nozero = amplitudes.nozero;
output.amplitude_nozero = amplitudes.amp_nozero;

% save output table
save([output_name '.mat'], 'output')
clear amplitudes

%% 6) unblind medication
% load and adjust the medication info
med_session = table2cell(readtable('treatments.xlsx'));
for a = 1 : size(med_session, 1)
    for b = 1 : length(session)
        if med_session{a, b+1}(1) == 'z'
            med_session{a, b+1} = 'placebo'; 
        else
            med_session{a, b+1} = 'alprazolam'; 
        end
    end
end
clear a b 

% change names of LW datasets
for p = 1:length(participant)
    for n = 1:length(session)
        this_session = med_session(participant(p), n+1);
        for t = 1:length(time)             
            for s = 1:length(stimulus) 
                name_old = [prefix_1 ' complete ' stimulus{s} ' YC' num2str(participant(p)) ' ' session{n} ' ' time{t} ' EMG'];
                
                % rename data with zero epochs
                load([name_old '.mat']); load([name_old '.lw6'], '-mat');                    
                name_new = [prefix_1 ' complete ' stimulus{s} ' YC' num2str(participant(p)) ' ' char(this_session) ' ' time{t} ' EMG'];
                header.name = name_new;
                save([name_new '.mat'], 'data');
                save([name_new '.lw6'], 'header');         
                
                % rename data without zero epochs
                load([prefix_2 ' ' name_old '.mat']); load([prefix_2 ' ' name_old '.lw6'], '-mat');                    
                name_new = [prefix_2 ' ' prefix_1 ' complete ' stimulus{s} ' YC' num2str(participant(p)) ' ' char(this_session) ' ' time{t} ' EMG'];
                header.name = name_new;
                save([name_new '.mat'], 'data');
                save([name_new '.lw6'], 'header'); 
            end
        end
    end
    message = ['Participant n.' num2str(participant(p)) ' finished.'];
    disp(message)
end
clear p n t s this_session name_old name_new data header

% change names in the output table
for p = 1:length(participant)
    for n = 1:length(session)
        this_session = med_session(participant(p), n+1);
        rows = (output.subject == participant(p) & categorical(output.session) == session{n});
        output.session(rows) = this_session;
    end
end
clear p n this_session rows

% sort and save the output table
output = sortrows(output, [1:4]);
output.medication = output.session; output.session = [];
output = output(:, [1, 12, 2:11]);
save([output_name '.mat'], 'output')

%% 7) plot group plots
% ----- spTMS MEPs - include zero epochs -----
% choose the data
data_visual = [];
for m = 1:length(medication)
    for t = 1:length(time)
        rows = (categorical(output.medication) == medication{m} & ...
            categorical(output.time) == time{t} & ...
            categorical(output.stimulus) == 'TS');
        data_i = output.amplitude_zero(rows);
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
figure_title = 'MEP amplitude - spTMS, 120 %rMT';
figure_name = 'MEP_TS';
yl = get(gca, 'ylim');
set(gca, 'xtick', 1:4, 'xticklabel', {'pre' 'post' 'pre' 'post'})
set(gca, 'Fontsize', 16)
title(figure_title, 'FontWeight', 'bold', 'FontSize', 18)
xlabel('time relative to medication'); ylabel('MEP (\muV)');
ylim([0, yl(2) + 200])
hold on

% add text
txt(1) = text(1.3, double(yl(2) + 100), 'placebo', 'fontsize', 16, 'color', colours2(2, :));
txt(2) = text(3.3, double(yl(2) + 100), 'alprazolam', 'fontsize', 16, 'color', colours2(4, :));
hold off

% save the figure       
savefig([figure_name '.fig'])
saveas(fig, [figure_name '.png'])

% update the counter
figure_counter = figure_counter + 1;     

% ----- changes in spTMS MEPs - include zero epochs -----
% choose the data
for p = 1:length(participant)
    data_visual_change(p, 1) = (data_visual(p, 2)/data_visual(p, 1)) * 100;
    data_visual_change(p, 2) = (data_visual(p, 4)/data_visual(p, 3)) * 100;
end

% plot group boxplot
col = colours2([2 4], :);
fig = figure(figure_counter);        
boxplot(data_visual_change, 'color', col)
hold on

% add no-change line
xl = get(gca, 'xlim');
h = line([xl(1) xl(2)], [100, 100], 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 2);
hold on

% plot the lines
for p = 1:length(participant)
    p_change(p) = plot([1 2], data_visual_change(p, [1 2]), '-o',...
        'Color', [0.75, 0.75, 0.75],...
        'MarkerSize', 10,...
        'MArkerEdge', 'none');
    hold on
end

% plot the markers
for b = 1:size(data_visual_change, 2)
    scat(b) = scatter(repelem(b, length(participant)), data_visual_change(:, b),...
        75, col(b, :), 'filled');
    hold on
end

% add parameters
figure_title = 'MEP amplitude change - spTMS, 120 %rMT';
figure_name = 'MEP_TS_change';
yl = get(gca, 'ylim');
set(gca, 'xtick', [1 2], 'xticklabel', {'placebo' 'alprazolam'})
set(gca, 'Fontsize', 16)
title(figure_title, 'FontWeight', 'bold', 'FontSize', 18)
xlabel('medication'); ylabel('MEP change (% baseline)');
ylim([yl(1), yl(2)])
hold off

% save the figure        
savefig([figure_name '.fig'])
saveas(fig, [figure_name '.png'])

% update the counter
figure_counter = figure_counter + 1;        

clear m t p b 
clear rows data_visual data_visual_change data_i p_placebo p_alprazolam p_change xl yl figure_title figure_name fig h scat txt col


               