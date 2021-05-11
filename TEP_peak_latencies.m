%% script to determine TOIs in TEPs of each single subject
% 1) imports baseline data of each participant and prepares the average
%       - crops the datasets according to preset params
%       - averages the data, creates a new header and saves for LW
% 2) plots data from selected electrodes and displays default TOI center
% 3) allows the user to manually identify the peak --> shifts the c to 
%    a new central value
% 4) saves the final TOI latencies of each subject into a table

%% parameters
clear all; clc;

% dataset
participant = 5;
stimulus = {'CS 80' 'TS 120' 'ppTMS 80120'};
medication = {'zyrtec' 'xanax'};
prefix = 'final_dataset avg avgchan';
prefix_new = 'all baseline avg crop';

% TOI
peaks = {'N15' 'P30' 'N45' 'P60' 'N100' 'P180'};
electrode = {'Cz', 'target'};
time_window = [-0.005, 0.4];
load('labels.mat')
load('peak_widths.mat')

% p = 5; s = 2; m = 1; 
%% core loop 
% prepare the outcome table
peak_latencies = table;
peak_latencies.subject = zeros(0); peak_latencies.channel = zeros(0); 
peak_latencies.peak = zeros(0); peak_latencies.latency = zeros(0); 

% start counting figures
figure_counter = 1;

% loop through participants
for p = 1:length(participant)
    
    % ------------------------- prepare the data ------------------------
    % import a random header
    load([prefix ' TS 120 YC' num2str(participant(1)) ' ' medication{1} ' pre M1.lw6'], '-mat')
    
    % loop through baseline datasets
    dataset_counter = 1;
    data_import = [];
    for s = 1:length(stimulus)
        for m = 1:length(medication)
            load([prefix ' ' stimulus{s} ' YC' num2str(participant(p)) ' ' medication{m} ' pre M1.mat'])
            % crop the data 
            data_import = cat(1, data_import, data(1, :, :, :, :, ...
                [((time_window(1) - header.xstart) / header.xstep) : ((time_window(2) - header.xstart) / header.xstep) + 1]));                       
            dataset_counter = dataset_counter + 1;
        end
    end
    
    % update the header
    header.name = [prefix_new ' YC' num2str(participant(p)) ' M1'];
    header.datasize(6) = length(squeeze(data_import));
    header.xstart = time_window(1);
    
    % average the data
    data = zeros(header.datasize);
    data(1, :, :, :, :, :) = mean(data_import, 1);
    
    % save new data + header
    save([header.name '.mat'], 'data');
    save([header.name '.lw6'], 'header');
    
    % -------------- identify the peaks for each electrode --------------
    for e = 1:length(electrode)
        % select the data from current electrode
        index = find(contains(labels, electrode(e))); 
        data_visual = double(squeeze(data(1, index, :, :, :, :)));
        
        % visualization parameters
        step = header.xstep;
        x = [time_window(1):step:time_window(2)];
        figure_name = ['YC' num2str(participant(p)) ' - average baseline TEPs - ' electrode{e}];
        
        % launch the figure
        axesHandles = [];
        fig = figure(figure_counter);
        axesHandles = [axesHandles subplot(4, 6, [7:24])];  
        plot(x, data_visual, 'b', 'LineWidth', 1)
        yl = ylim;
        xlim([time_window(1), time_window(2)])
        hold on
        
        rectangle('Position', [0, yl(1), 0.01, yl(2)-yl(1)], 'FaceColor', [0.9020    0.9020    0.9020], 'EdgeColor', 'none')
        title(figure_name, 'FontSize', 16)
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        plot(x, data_visual, 'b', 'LineWidth', 3)
        line([0, 0], get(gca,'ylim'), 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2)
        line(x, zeros(length(x)), 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 1)
        hold on
        
        % loop throught he peaks and set all individual TOIs
        for peak = 1:length(peaks)
            % identify the default peak latency and span
            switch electrode{e}
                case 'Cz'
                        center = peak_widths.latency(peak);
                        span = peak_widths.width(peak)/1000;                  
                case 'target'
                        center = peak_widths.latency(peak + length(peaks));
                        span = peak_widths.width(peak + length(peaks))/1000;   
            end
            
            % plot the default peak center 
            subplot(4, 6, [7:24])
            line([center, center], get(gca,'ylim'), 'Color', 'red', 'LineWidth', 1, 'LineStyle', ':');    
            hold on
            
            % add the topoplot 
            map_lims = get(axesHandles(1), 'ylim');
            axesHandles = [axesHandles subplot(4, 6, peak)];
            TEP_topoplot(header, data, center, map_lims);
            title(peaks{peak}, 'FontSize', 16)
            hold on
            
            % let the topoplot change with cursor x position
            set(fig, 'WindowButtonMotionFcn', {@mouse_move2, axesHandles, header, data, peak});  
            
            % wait for a click
            pos_x = get_position(axesHandles);
            
            % update the figure
            set (fig, 'WindowButtonMotionFcn', '');
            subplot(4, 6, [7:24])
            line([pos_x, pos_x], get(gca,'ylim'), 'Color', 'red', 'LineWidth', 2);
            subplot(4, 6, peak) 
            cla(axesHandles(peak+1))
            TEP_topoplot(header, data, pos_x, map_lims);
            hold on
            
            % save current center value to the table
            latency_i = table;
            latency_i.subject = participant(p); latency_i.channel = electrode(e); 
            latency_i.peak = peaks(peak); latency_i.latency = pos_x; 
            
            peak_latencies = [peak_latencies; latency_i];            
        end
        
        % save the figure
        savefig([figure_name '.fig'])
        
        % update the figure counter
        figure_counter = figure_counter + 1;
        
    end
end

% save the latencies
save('peak_latencies.mat', 'peak_latencies')
