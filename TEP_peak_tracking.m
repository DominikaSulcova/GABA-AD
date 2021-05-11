%% parameters
clear all; clc;

% dataset
participant = 6:20;
stimulus = {'CS 80'};  
medication = {'zyrtec' 'xanax'};
time = {'pre'};
prefix = 'bl final_dataset avg avgchan';
prefix_new = 'tracked_Cz';

% peaks
peaks = {'N15' 'P30' 'N45' 'P60' 'N100' 'P180'};
electrode = {'Cz'};
load('labels.mat')
load('peak_widths.mat')
buffer_percent = 0.5;       % a margin of the window for peak visualisation
% k = 2; p = 5; s = 1; m = 1; t = 1;

% figure numbering 
figure_counter = 1;

%% load the data and identify the peak in Cz
% check for the table of latencies
if exist('CS_latencies.mat') == 0
    CS_latencies = table;
    CS_latencies.subject = zeros(0); CS_latencies.stimulus = zeros(0); CS_latencies.medication = zeros(0); CS_latencies.time = zeros(0);
    CS_latencies.peak = zeros(0); CS_latencies.latency = zeros(0); 
else
    load('CS_latencies.mat')
end

% loop through the participants
for k = 2:length(peaks)
    for p = 1:length(participant)
        for m = 1:length(medication)
            for s = 1:length(stimulus)
                for t = 1:length(time)                
                    % load the da data and the header
                    load([prefix ' ' char(stimulus(s)) ' YC' num2str(participant(p)) ' ' char(medication(m)) ' ' char(time(t)) ' M1.lw6'],'-mat');      
                    load([prefix ' ' char(stimulus(s)) ' YC' num2str(participant(p)) ' ' char(medication(m)) ' ' char(time(t)) ' M1.mat']);

                    % deal with filenames
                    switch stimulus{s}
                        case 'CS 80'
                            stim = 'CS';
                        case 'TS 120'
                            stim = 'TS';
                    end                
                    switch medication{m}
                        case 'xanax'
                            med = 'alprazolam';
                        case 'zyrtec'
                            med = 'placebo';
                    end
                    
                    % figure params                
                    figure_name = ['YC' num2str(participant(p)) ' - ' med ', ' peaks{k}];
                    figure_center = peak_widths.latency(k);
                    span = ((1 + buffer_percent) * peak_widths.width(k)) / 1000;
                    
                    finish = 0;
                    while finish == 0
                        % identify the TOI of current peak                    
                        x1 = ceil((figure_center - span/2 - header.xstart) / header.xstep);
                        x2 = ceil((figure_center + span/2 - header.xstart) / header.xstep);
                        x = (figure_center - span/2) : header.xstep : (figure_center + span/2);

                        % prepare data and header for visualization
                        data_visual = data(1, :, :, :, :, [x1 : x2]);
                        header_visual = header;
                        header_visual.datasize(6) = length(data_visual);  
                        header_visual.xstart = figure_center - span/2;

                        % launcg the figure                    
                        axesHandles = [];
                        fig = figure(figure_counter);   
                        axesHandles = [axesHandles subplot(3, 3, [4:9])];  
                        plot(x, double(squeeze(data_visual(1, find(contains(labels, electrode{1})), :, :, :, :))), 'LineWidth', 2)
                        xlim([(figure_center - span/2), (figure_center + span/2)])
                        title(figure_name, 'FontSize', 16)
                        hold on                

                        % plot the line at the center
                        h = line([figure_center, figure_center], get(gca,'ylim'), 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--'); 
                        hold on    

                        % plot the central topography 
                        map_lims = get(axesHandles(1), 'ylim');
                        axesHandles = [axesHandles subplot(3, 3, 2)];
                        TEP_topoplot(header_visual, data_visual, figure_center, map_lims);
                        hold on

                        % make the topography change with mouse movement 
                        set(fig, 'WindowButtonMotionFcn', {@mouse_move, axesHandles, header_visual, data_visual});              

                        % choose the peak position
                        pos_x = get_position(axesHandles);             

                        % update the figure
                        set (fig, 'WindowButtonMotionFcn', '');
                        subplot(3, 3, [4:9])
                        set(h, 'XData', [pos_x, pos_x], 'LineStyle', '-');
                        subplot(3, 3, 2) 
                        cla(axesHandles(2))
                        TEP_topoplot(header_visual, data_visual, pos_x, map_lims);
                        hold on

                        pause(2)

                        % ask for approval
                        answer = questdlg('Do you want to proceed?', figure_name,...
                            'Yes', 'No, I want to keep fiddling', 'Yes');

                        % switch action
                        switch answer
                            case 'Yes'
                                % close the figure
                                close(fig)

                                % exit the while loop
                                finish = 1;

                            case 'No, I want to keep fiddling'
                                % assign previous center
                                figure_center = pos_x;  
                                
                                % close the figure
                                close(fig)                                
                        end
                    end                    

                    % index the peak timepoint  
                    n_peak = ceil((pos_x - header.xstart) / header.xstep); 

                    % loop through the data and subtract the value at the peak latency from all
                    % timepoints
                    sub_data = data;
                    topo_vector = squeeze(data(1, :, 1, 1, 1, n_peak));                                             
                    for i = 1:header.datasize(2)    % all channels
                        for n = 1:length(data)      % all datapoints 
                            sub_data(1, i, :, :, :, n) = data(1, i, :, :, :, n) - topo_vector(i);
                        end
                    end

                    % crop the data
                    x_out1 = ceil((pos_x - span/2 - header.xstart) / header.xstep);
                    x_out2 = ceil((pos_x + span/2 - header.xstart) / header.xstep);
                    sub_data = sub_data(1, :, :, :, :, [x_out1 : x_out2]);
                    data = data(1, :, :, :, :, [x_out1 : x_out2]);

                    % create a header
                    header.datasize(6) = length(data);
                    header.xstart = - span/2;           % put the peak latency to 0
                    header.events.code = 'peak_Cz';     % rename the event

                    % save the original cropped data
                    name = [prefix_new ' ' peaks{k} ' original ' char(stimulus(s)) ' YC' num2str(participant(p)) ' ' char(medication(m)) ' ' char(time(t)) ' M1'];
                    header.name = name; 
                    save([name '.mat'], 'data');
                    save([name '.lw6'], 'header');

                    % save the subtracted data
                    name = [prefix_new ' ' peaks{k} ' subtracted ' char(stimulus(s)) ' YC' num2str(participant(p)) ' ' char(medication(m)) ' ' char(time(t)) ' M1'];
                    header.name = name;
                    data = sub_data;
                    save([name '.mat'], 'data');
                    save([name '.lw6'], 'header'); 

                    % fill in the table of latencies
                    latency_i = table;
                    latency_i.subject = participant(p);
                    latency_i.stimulus = {stim};
                    latency_i.medication = {med};
                    latency_i.time = time(t);
                    latency_i.peak = peaks(k);
                    latency_i.latency = pos_x; 

                    CS_latencies = [CS_latencies; latency_i];

                    % save latencies
                    save('CS_latencies.mat', 'CS_latencies')

                    % update the counter                                
                    figure_counter = figure_counter + 1;             
                end
            end
        end
    end
end
