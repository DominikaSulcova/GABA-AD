%% CALCULATE MEAN AMPLITUDE OF TEP PEAKS
% Written by Dominika for GABA-AD project (2020)
% 
% 1) Uses TOI widths calculated from baseline grand average TEPs (Cz) and centers
%    them around individual average baseline peak latencies (selected manually
%     - see script GABA_TEP_latencies.m)
% 2) Calculates mean of <percent> datapoints of each TOI
% 3) Indentifies central latency of these averaged datapoints
% 4) Plots TEPs of selected electrodes with the visualisation of AOCs and
%    appropriate topoplots; saves the figure
% 5) Saves extracted amplitudes and latencies in a common table 'GABA_TEP_<subject group>'

%% parameters
clear all; clc;

% dataset
subj_group = 'YC';
participant = 1:20;
stimulus = {'TS' 'ppTMS'};
medication = {'placebo' 'alprazolam'};
time = {'pre' 'post'};
prefix = 'eois crop avg final_dataset';
load('labels.mat')

% TOI
electrode = {'eoi_P30' 'eoi_N45_central' 'eoi_P60' 'eoi_N100' 'eoi_P180'};
peaks = {'P30' 'N45' 'P60' 'N100' 'P180'};
time_window = [-0.005, 0.3];
load('peak_widths.mat'); 
PW = peak_widths(13:height(peak_widths), :);

% amplitude
percent = 25;

% p = 5; s = 2; m = 1; t = 1; e = 1;

%% extract mean amplitudes from current subjects
% prepare the main table of TEP results
outcome_name = ['GABA_TEP_' subj_group];
if exist([outcome_name '.mat']) == 0
    TEPs = table;
    TEPs.subject = zeros(0); TEPs.stimulus = zeros(0); TEPs.medication = zeros(0); TEPs.time = zeros(0);
    TEPs.channel = zeros(0); TEPs.peak = zeros(0); TEPs.amplitude = zeros(0); TEPs.latency = zeros(0); 
else
    load('TEPs_final.mat')
end

% set the figure counter
figure_counter = 1; 

% loop through all subjects
n = length(electrode);
for p = 1:length(participant) 
    for m = 1:length(medication)  
        for s = 1:length(stimulus)
            for t = 1:length(time)
                % prepare outcome table
                TEPs_i = table;
                
                % setup names
                figure_name = ['YC' num2str(participant(p)) ': ' stimulus{s} ', ' medication{m} ', ' time{t} ' medication'];
                file_name = ['YC' num2str(participant(p)) '_' stimulus{s} '_' medication{m} '_' time{t}];
                
                % loop through the electrodes
                for e = 1:length(electrode)  
                    % in case of TS pre medication:
                    if s == 1 & t == 1 & e == 1
                        % use default PW
                        PW_i = PW;
                        PW_i.electrode{2} = 'eoi_N45_central';
                        
                        % alert the operator at the beginning of each
                        % baseline TS
                        if  e == 1                       
                            load train
                            sound(y,Fs)
                        end
                    end
                    
                    % load the da data and the header
                    load([prefix ' ' char(stimulus(s)) ' YC' num2str(participant(p)) ' ' char(medication(m)) ' ' char(time(t)) ' M1.lw6'],'-mat');     
                    load([prefix ' ' char(stimulus(s)) ' YC' num2str(participant(p)) ' ' char(medication(m)) ' ' char(time(t)) ' M1.mat']);                
                                        
                    % set up visualisation parameters
                    x = [time_window(1):header.xstep:time_window(2)];
                    x_start = (time_window(1)-header.xstart)/header.xstep;
                    x_end = (time_window(2)-header.xstart)/header.xstep;
                    col = [0,0.451,0.740];
                    
                    % prepare data for visualization
                    position = find(contains(labels, electrode{e})); 
                    data_visual = double(squeeze(data(1, position, :, :, :, x_start:x_end))');

                    % define default TOI
                    if e == length(electrode) & strcmp(stimulus{s}, 'CS')
                        rows = (categorical(PW_i.electrode) == 'eoi_P180');
                    else
                        rows = (categorical(PW_i.electrode) == electrode{e});
                    end
                    center = PW_i.latency(rows);
                    span = PW_i.width(rows)*header.xstep;
                    
                    % set the final TOI
                    finish = 0;
                    while finish == 0;
                        % launch the figure
                        fig = figure(1);  
                        axesHandles = [];
                        
                        % plot the data timeseries
                        axesHandles = [axesHandles subplot(4, 6, [7:24])]; 
                        plot(x, data_visual, 'b', 'LineWidth', 1)
                        yl = ylim;
                        xlim([time_window(1), time_window(2)])
                        rectangle('Position', [0, yl(1), 0.01, yl(2)-yl(1)], 'FaceColor', [0.9020    0.9020    0.9020], 'EdgeColor', 'none')
                        title(peaks{e}, 'FontSize', 16)
                        set(gcf,'units','normalized','outerposition',[0 0 1 1])
                        hold on

                        % set limits for topoplot colourmap
                        map_lims = get(axesHandles(1), 'ylim');
                        
                        % visualize default peak TOI
                        subplot(4, 6, [7:24])
                        rectangle('Position', [center - span/2, yl(1), span, yl(2)-yl(1)], 'FaceColor', [1,0.7608,0.7608], 'EdgeColor', 'none')

                        % define the polarity of current peak 
                        if mod(e, 2) == 1
                            polarity = 'positive';
                        else
                            polarity = 'negative';
                        end  
                        
                        % calculate mean amplitude
                        [mean_amplitude, averaged_x, averaged_data] = GABA_TEP_meanamp(data_visual, polarity, center, span, percent, header.xstep, time_window(1)); 

                        % calculate  real peak latency (round up)
                        central_latency = averaged_x(ceil(length(averaged_x)/2));

                        % update the figure
                        subplot(4, 6, [7:24])
                        for a = 1:length(averaged_x)
                            line([averaged_x(a), averaged_x(a)], [0, averaged_data(a)-0.05], 'Color', 'red', 'LineWidth', 1)
                            hold on
                        end
                        
                        % add the topoplot    
                        axesHandles = [axesHandles subplot(4, 6, 1)];
                        GABA_TEP_topoplot(header, data, central_latency, map_lims);
                        hold on 
                        
                        % replot the data to make it visible
                        subplot(4, 6, [7:24])
                        plot(x, data_visual, 'b', 'LineWidth', 2.5)
                        line([0, 0], get(gca,'ylim'), 'LineStyle', '--', 'Color', col, 'LineWidth', 2.5)
                        line(x, zeros(length(x)), 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 1)
                        hold on
                        
                        % ask for approval
                        if s == 1 & t == 1 
                            answer = questdlg('Do you want to proceed?', figure_name,...
                                'Yes, extract outcome values.', 'No, I will adjust TOI manually.', 'Yes, extract outcome values.');
                        else
                            answer = 'Yes, extract outcome values.';
                        end
                        
                        % switch action
                        switch answer
                            case 'Yes, extract outcome values.'
                                % close the figure
                                close(fig)
                                
                                % exit the while loop
                                finish = 1;

                            case 'No, I will adjust TOI manually.'
                                % assign previous center and span
                                choose_center = center;  
                                choose_span = 2 * span;  

                                % identify the limits for visualisation of current peak
                                choose_x1 = ceil((choose_center - choose_span/2 - header.xstart) / header.xstep);
                                choose_x2 = ceil((choose_center + choose_span/2 - header.xstart) / header.xstep);
                                choose_x = (choose_center - choose_span/2) : header.xstep : (choose_center + choose_span/2);

                                % prepare data and header for visualization
                                choose_data = double(squeeze(data(1, position, :, :, :, [choose_x1 : choose_x2])));
                                choose_header = header;
                                choose_header.datasize(6) = length(choose_data);  
                                choose_header.xstart = choose_center - choose_span/2;
                                
                                % check if vector size matches
                                if length(choose_data) ~= length(choose_x)
                                    diff = length(choose_data) - length(choose_x);
                                    if diff > 0
                                        choose_data = choose_data(1:end - diff);
                                    elseif diff < 0
                                        choose_x = choose_x(1:end + diff);
                                    end
                                end

                                % launch the choosing figure                 
                                choose_figure_name = ['Choose manually peak ' peaks{e}];
                                choose_axesHandles = [];
                                choose_fig = figure(2);   
                                choose_axesHandles = [choose_axesHandles subplot(3, 3, [4:9])];  
                                plot(choose_x, choose_data, 'LineWidth', 2)
                                xlim([(choose_center - choose_span/2), (choose_center + choose_span/2)])
                                title(choose_figure_name, 'FontSize', 16)
                                hold on                

                                % plot the line at the center
                                h = line([choose_center, choose_center], get(gca,'ylim'), 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--'); 
                                hold on    

                                % plot the central topography 
                                choose_map_lims = get(choose_axesHandles(1), 'ylim');
                                choose_axesHandles = [choose_axesHandles subplot(3, 3, 2)];
                                GABA_TEP_topoplot(header, data, choose_center, choose_map_lims);
                                hold on

                                % make the topography change with mouse movement 
                                set(choose_fig, 'WindowButtonMotionFcn', {@mouse_move, choose_axesHandles, header, data,});              

                                % choose the peak position
                                pos_x = get_position(choose_axesHandles);  
                                set(choose_fig, 'WindowButtonMotionFcn', '')

                                % update the figure
                                set (choose_fig, 'WindowButtonMotionFcn', '');
                                subplot(3, 3, [4:9])
                                set(h, 'XData', [pos_x, pos_x], 'LineStyle', '-');
                                subplot(3, 3, 2) 
                                cla(choose_axesHandles(2))
                                GABA_TEP_topoplot(header, data, pos_x, map_lims);
                                hold off

                                % update the central latency
                                center = pos_x;

                                % close the choosing figure
                                pause(5)
                                close(choose_fig)

                                % close the the main figure
                                close(fig)
                        end
                    end
                    
                    % update PW_i for future datasets
                    if s == 1 & t == 1 
                        PW_i.latency(rows) = center;
                    end
                    
                    % record outcome variables
                    TEPs_i.subject(e, 1) = participant(p);
                    TEPs_i.stimulus(e, 1) = stimulus(s);
                    TEPs_i.medication(e, 1) = medication(m);
                    TEPs_i.time(e, 1) = time(t);
                    TEPs_i.channel(e, 1) = electrode(e);
                    TEPs_i.peak(e, 1) = peaks(e);
                    TEPs_i.amplitude(e, 1) = mean_amplitude; 
                    TEPs_i.latency(e, 1) = central_latency;  
                    
                    % launch summary figure 
                    fig = figure(figure_counter);  
                    axesHandles = [];
                    
                    % plot the data timeseries
                    b = (1 + 3 * (e - 1));
                    c = (2 + 3 * (e - 1));  
                    d = (3 + 3 * (e - 1));  
                    axesHandles = [axesHandles subplot(3, 6, [b, c])];  
                    plot(x, data_visual, 'b', 'LineWidth', 1)
                    yl = ylim;
                    xlim([time_window(1), time_window(2)])
                    rectangle('Position', [0, yl(1), 0.01, yl(2)-yl(1)], 'FaceColor', [0.9020    0.9020    0.9020], 'EdgeColor', 'none')
                    title(peaks{e}, 'FontSize', 16)
                    set(gcf,'units','normalized','outerposition',[0 0 1 1])
                    hold on
                    
                    % set limits for topoplot colourmap
                    map_lims = get(axesHandles(1), 'ylim');

                    % visualize final peak TOI
                    subplot(3, 6, [b, c])
                    rectangle('Position', [center - span/2, yl(1), span, yl(2)-yl(1)], 'FaceColor', [1,0.7608,0.7608], 'EdgeColor', 'none')
                    
                    % update the figure
                    subplot(3, 6, [b, c])
                    for a = 1:length(averaged_x)
                        line([averaged_x(a), averaged_x(a)], [0, averaged_data(a)-0.05], 'Color', 'red', 'LineWidth', 1)
                        hold on
                    end

                    % add the topoplot    
                    axesHandles = [axesHandles subplot(3, 6, d)];
                    GABA_TEP_topoplot(header, data, central_latency, map_lims);
                    hold on 

                    % replot the data to make it visible
                    subplot(3, 6, [b, c])
                    plot(x, data_visual, 'b', 'LineWidth', 2, 'Color', col)
                    line([0, 0], get(gca,'ylim'), 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2)
                    line(x, zeros(length(x)), 'LineStyle', ':', 'Color', [0, 0, 0], 'LineWidth', 1)
                    suptitle(figure_name)
                    hold off                    
                end
                
                % append the results to the main table, save
                TEPs = [TEPs; TEPs_i];
                save([outcome_name '.mat'], 'TEPs')
                
                % save the summary figure, close
                savefig(fig, [file_name '.fig'])
                close(fig)
                
                % update the figure counter
                figure_counter = figure_counter + 1;                      
                    
            end
        end
    end
    
    % play a celebratory sound at the end of this friggin participant
    tune = load('handel.mat');
    sound(tune.y, tune.Fs)
end



