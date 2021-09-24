%% GABA-AD: TMS-EVOKED POTENTIALS
% Written by Dominika for GABA-AD project (2020-21)
% 
% % Colection of scripts to analyse preprocessed TMS-EEG data: 
%       * all data and figures are saved in a new folder *
% 1) load the data
%       - loads individual averaged data and trims it in a predefined time
%       window
%       - saves the data from all subjects and conditions in a 6D matrix:
%       2 (medication) X 2 (time) X number of stimuli  X number of subjects X
%       number of channels X number of timepoints
% 
% 2) preliminary visualization of TEPs 
%       - plots baseline TEPs from selected electrode, both conditions in
%       one plot --> mean + CI
% 
% 3) calculate global mean field power from grand average
%       - averages baseline and post-medication data from both sessions
%       - calculates overall GMFP and plots it, adds topoplots for peak
%       times, saves the figure
%       - automatically identifies local peaks within [0.01 0.25]s and
%       saves peak latencies in the outcome structure
% 
% 4) define peak widths using grand average GMFP
%       - for each stimulus separately, plots and saves the figure 
%       - uses function findpeaks awith the 'halfheight' option
%       - saves peak widths in the outcome structure
% 
% 5) determine electrodes of interest
%       - uses GMFP peak latencies as default search points and signal from Cz 
%       as starting peak point to track peak progress
%       - averages across subjects
%       - identifies 3 electrodes with the highest average amplitude for
%       each peak, plots the  




%% parameters
clear all; clc

% ----- adjustable parameters -----
% dataset
group = 'YC';
target = 'AG';
participant = 1:20;
medication = {'placebo' 'alprazolam'}; 
time = {'pre' 'post'};

% visualization 
time_window = [-0.05, 0.3];
z = 1.96;
alpha = 0.2;
% --------------------------------

% dataset parameters
prefix = 'eoi avg bl icfilt ica26 crop but fft-notchfilt prefilt prea28 visual reref ds art-sup bl ep dc chan-select';
load([prefix ' ' group num2str(participant(1)) ' ' medication{1} ' ' target ' ' time{1} '.lw6'], '-mat')
labels.AG = {header.chanlocs.labels};

% visualization calculated params
figure_counter = 1;
xstep = header.xstep; 
x = [time_window(1):xstep:time_window(2)];
x_start = (time_window(1) - header.xstart)/xstep;
x_end = (time_window(2) - header.xstart)/xstep;
load('colours.mat'); 

% create output folders
filename = ['GABA_' group '_' target];
folderpath = [pwd '\' filename];
if ~exist(folderpath) 
    mkdir(folderpath)
    mkdir(folderpath, [filename '_figures'])
end     

%% 1) extract individual data
% load data based on the prefix + conditions
n = length(labels.TS) - length(labels.CS);
% m = 1;t = 1; s = 1; p = 1;
for m = 1:length(medication)
    for t = 1:length(time)
        for p = 1:length(participant)
            % load individual dataset
            load([prefix ' ' group num2str(participant(p)) ' ' medication{m} ' ' target ' ' time{t} '.mat'])

            % append the data in the data matrix
            GABA_data(m, t, p, :, :) = squeeze(data(:, :, :, :, :, x_start:x_end));                        
        end
    end
end

clear m t s p data n
disp(['Data import finished. Datasize: ' num2str(size(GABA_data))])

% save dataset to the global MATLAB file
if ~exist([folderpath '\' filename '.mat']) 
    save([folderpath '\' filename '.mat'], 'GABA_data');
else
    save([folderpath '\' filename '.mat'], 'GABA_data', '-append');
end

%% 2) preliminary visualization 
% ----- decide output parameters -----
electrode = {'target'};
% ------------------------------------

% prepare average statistics
for m = 1:length(medication)
    for t = 1:length(time)
        for e = 1:size(GABA_data, 4)
            for i = 1:size(GABA_data, 5)
                data_mean(m, t, e, i) = mean(squeeze(GABA_data(m, t, :, e, i)));
                data_CI(m, t, e, i) = (std(squeeze(GABA_data(m, t, :, e, i)))/sqrt(length(participant))) * z;
            end
        end
    end
end
clear m t e i

% plot the baseline data both sessions
for e = 1:length(electrode) 
    % identify the electrode
    e_n = find(contains(labels, electrode{e}));
    
    % prepare data
    data_visual = squeeze(data_mean(:, 1, e_n, :));
    CI_visual = squeeze(data_CI(:, 1, e_n, :));

    % launch the figure
    fig = figure(figure_counter);
    hold on

    % set limits of the figure
    plot(x, data_visual(1, :) + CI_visual(1, :), 'b:', 'LineWidth', 0.5)
    plot(x, data_visual(1, :) - CI_visual(1, :), 'b:', 'LineWidth', 0.5)
    yl = get(gca, 'ylim'); yl(1) = yl(1) - 0.2; yl(2) = yl(2) + 0.3;
    clf, hold on

    % shade interpolated interval 
    plot(x, data_visual(1, :), 'b:', 'LineWidth', 0.5);
    rectangle('Position', [0, yl(1), 0.01, yl(2) - yl(1)], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')

    % loop through datasets to plot
    for m = 1:length(medication)        
        P(m) = plot(x, data_visual(m, :), 'Color', colours((m - 1)*2 + 1, :), 'LineWidth', 2.5);
        F(m) = fill([x fliplr(x)],[data_visual(m, :) + CI_visual(m, :) fliplr(data_visual(m, :) - CI_visual(m, :))], ...
            colours((m - 1)*2 + 1, :), 'FaceAlpha', alpha, 'linestyle', 'none');
    end

    % add other parameters
    title([target ', ' electrode{e} ' electrode: baseline TEP'])
    xlabel('time (s)')
    ylabel('amplitude (\muV)')
    set(gca, 'FontSize', 14)
    xlim(time_window)
    ylim(yl)
    line([0, 0], get(gca,'ylim'), 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)
    lgd = legend(P, medication, 'Location', 'southeast');
    lgd.FontSize = 14;
    hold off

    % name and save figure
    figure_name = ['TEP_' target '_bl_' electrode{e}];
    savefig([folderpath '\' filename '_figures\' figure_name '.fig'])
    saveas(fig, [folderpath '\' filename '_figures\' figure_name '.png'])

    % update figure counter
    figure_counter = figure_counter + 1 ;
end
clear m e e_n fig lgd data_visual CI_visual figure_name P F yl

% plot baseline vs. post-medication 
for e = 1:length(electrode) 
    for m = 1:length(medication)
        % identify the electrode
        e_n = find(contains(labels, electrode{e}));

        % prepare data
        data_visual = squeeze(data_mean(m, :, e_n, :));
        CI_visual = squeeze(data_CI(m, :, e_n, :));

        % launch the figure
        fig = figure(figure_counter);
        hold on

        % set limits of the figure
        plot(x, data_visual(1, :) + CI_visual(1, :), 'b:', 'LineWidth', 0.5)
        plot(x, data_visual(1, :) - CI_visual(1, :), 'b:', 'LineWidth', 0.5)
        yl = get(gca, 'ylim'); yl(1) = yl(1) - 0.2; yl(2) = yl(2) + 0.3;
        clf, hold on

        % shade interpolated interval 
        plot(x, data_visual(1, :), 'b:', 'LineWidth', 0.5);
        rectangle('Position', [0, yl(1), 0.01, yl(2) - yl(1)], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')

        % loop through datasets to plot
        for t = 1:length(time)        
            P(t) = plot(x, data_visual(t, :), 'Color', colours((m - 1)*2 + t, :), 'LineWidth', 2.5);
            F(t) = fill([x fliplr(x)],[data_visual(t, :) + CI_visual(t, :) fliplr(data_visual(t, :) - CI_visual(t, :))], ...
                colours((m - 1)*2 + t, :), 'FaceAlpha', alpha, 'linestyle', 'none');
        end

        % add other parameters
        title([target ', ' electrode{e} ' electrode: ' medication{m}])
        xlabel('time (s)')
        ylabel('amplitude (\muV)')
        set(gca, 'FontSize', 14)
        xlim(time_window)
        ylim(yl)
        line([0, 0], get(gca,'ylim'), 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)
        lgd = legend(P, {'baseline' 'post medication'}, 'Location', 'southeast');
        lgd.FontSize = 14;
        hold off

        % name and save figure
        figure_name = ['TEP_' target '_' medication{m} '_' electrode{e}];
        savefig([folderpath '\' filename '_figures\' figure_name '.fig'])
        saveas(fig, [folderpath '\' filename '_figures\' figure_name '.png'])

        % update figure counter
        figure_counter = figure_counter + 1 ;
    end
end
clear m e t e_n fig lgd data_visual CI_visual figure_name P F yl

% save dataset to the global MATLAB file
save([folderpath '\' filename '.mat'], 'data_mean', 'data_CI', '-append');
clear electrode 

%% 3) GMFP - peak identification
% ----- decide output parameters -----
labeled = 'off';
max_peaks = 6;
% ------------------------------------

% calculate global GMFP (exclude target and eoi channels)
for m = 1:length(medication)
    for t = 1:length(time)
        GABA_GMFP(m, t, :) = std(squeeze(data_mean(m, t, 1:30, :)), 1);  
    end 
end
clear m t 

% remove excessive channels in the header
header.chanlocs(31:end) = [];
header.datasize(2) = 30;  

% pool all conditions together
for i = 1:size(GABA_GMFP, 3)
    data_visual(i) = mean(GABA_GMFP(:, :, i), 'all');
end

% launch the figure
fig = figure(figure_counter);
hold on

% extract peak latencies
h_axis(1) = subplot(3, max_peaks, [1 : 2*max_peaks]);
GABA_peaks = gmfp_plot(x, data_visual, time_window, xstep, labeled, 'max_peaks', max_peaks);
title(['M1 - ' stimulus{s} ', all conditions: GMFP'], 'fontsize', 16, 'fontweight', 'bold')

% choose data for topoplots 
for e = 1:30
    for i = 1:size(data_mean, 5)
        data_topoplot(1, e, 1, 1, 1, i) = mean(data_mean(:, :, e, i), 'all');
    end
end

% add topoplots
for t = 1:length(GABA_peaks)
    % plot the topoplot
    h_axis(1 + t) = subplot(3, max_peaks, 2*max_peaks + t);
    topo_plot(header, data_topoplot, GABA_peaks(t), time_window(1), [-2, 2])

    % shift down
    pos = get(h_axis(1 + t), 'Position');
    pos(2) = pos(2) - 0.05;
    set(h_axis(1 + t), 'Position', pos);

    % add timing
    text(-0.3, -0.8, sprintf('%1.0f ms', GABA_peaks(t)*1000), 'Color', [1 0 0], 'FontSize', 14)
end
hold off

% save figure
figure_name = ['GMFP_' target];
savefig([folderpath '\' filename '_figures\' figure_name '.fig'])
saveas(fig, [folderpath '\' filename '_figures\' figure_name '.png'])

% update figure counteer
figure_counter = figure_counter + 1;
clear i t e data_visual data_topoplot fig figure_name pos h_axis GABA_peaks

% append new variables to the general MATLAB file
save([folderpath '\' filename '.mat'], 'GABA_GMFP', '-append');
clear labeled max_peaks 

%% 4) peak widths                    
%% 5) EOIs
%% 6) amplitude extraction
%% functions
function peak_x = gmfp_plot(x, y, time_window, xstep, labeled, varargin)
% check whether to plot labels (default)
if ~isempty(varargin)
    a = find(strcmpi(varargin, 'max_peaks'));
    if ~isempty(a)
        max_peaks = varargin{a + 1};
    end
end

% launch the figure  
plot(x, y)
yl = get(gca, 'ylim');
cla

% plot interpolated part
hold on
xlim(time_window)
rectangle('Position', [0, yl(1), 0.01, yl(2) - yl(1)], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')

% plot data, mark TMS stimulus
plot(x, y, 'Color', [0 0 0], 'LineWidth', 2.5)
line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)

% find peaks 
[pks, locs] = findpeaks(y, 'MinPeakDistance', 10, 'MinPeakProminence', 0.015);
for a = 1:length(locs)
    if time_window(1) + locs(a)*xstep <= 0.015
        idx(a) = false;
    elseif time_window(1) + locs(a)*xstep > 0.220
        idx(a) = false;        
    else
        idx(a) = true;
    end
end
pks = pks(idx); locs = locs(idx);
if length(pks) > max_peaks
    pks = pks(1:max_peaks); 
    locs = locs(1:max_peaks);
end

% calculate peak coordinations
for a = 1:length(locs)
    peak_x(a) = time_window(1) + locs(a)*xstep;
    peak_y(a) = pks(a);
end
peak_y = double(peak_y);

% plot peaks
for a = 1:length(locs)
    plot(peak_x(a), peak_y(a), ...
        'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none');
    line([peak_x(a), peak_x(a)], [yl(1), peak_y(a)], 'LineStyle', ':', 'Color', [1, 0, 0], 'LineWidth', 1.5)
    
    % label with latency (ms)
    if strcmp(labeled, 'on') 
        text(peak_x(a), peak_y(a) + 0.15, sprintf('%1.0fms', peak_x(a)*1000), 'Color', [1 0 0], 'FontSize', 14)
    end
end

% add parameters
set(gca, 'fontsize', 14)
ylim(yl)
xlabel('time (s)')
ylabel('potential (\muV)')
end
function [peak_x, peak_y] = gmfp_peaks(y, time_window, xstep, varargin)
% check whether to plot labels (default)
if ~isempty(varargin)
    a = find(strcmpi(varargin, 'max_peaks'));
    if ~isempty(a)
        max_peaks = varargin{a + 1};
    end
end

% find peaks 
[pks, locs] = findpeaks(y, 'MinPeakDistance', 10, 'MinPeakProminence', 0.015);
for a = 1:length(locs)
    if time_window(1) + locs(a)*xstep <= 0.015
        idx(a) = false;
    elseif time_window(1) + locs(a)*xstep > 0.220
        idx(a) = false;        
    else
        idx(a) = true;
    end
end
pks = pks(idx); locs = locs(idx);
if length(pks) > max_peaks
    pks = pks(1:max_peaks); 
    locs = locs(1:max_peaks);
end

% calculate peak coordinations
for a = 1:length(locs)
    peak_x(a) = time_window(1) + locs(a)*xstep;
    peak_y(a) = pks(a);
end
peak_y = double(peak_y);
end
function [p1, p2, f] = gmfp_plot_diff(x, y, time_window, colours)
% calculate the time difference
for g = 1:size(y, 2)
    y(3, g) = y(2, g) - y(1, g);
end

% launch the figure 
hold on
plot(x, y(1, :))
plot(x, y(3, :))
yl = get(gca, 'ylim'); yl = [yl(1) - 0.1 , yl(2) + 1];
clf

% plot interpolated part
hold on
xlim(time_window)
rectangle('Position', [0, yl(1), 0.01, yl(2) - yl(1)], 'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none')

% plot data
p1 = plot(x, y(1, :), 'Color', [0 0 0], 'LineWidth', 2.5);
p2 = plot(x, y(2, :), 'Color', colours(2, :) , 'LineWidth', 2.5);
f = fill([x fliplr(x)],[y(3, :) fliplr(zeros(1, length(x)))], colours(1, :), 'linestyle', 'none');
    
% add lines    
line([0, 0], yl, 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5)
line(time_window, [0, 0], 'Color', [0, 0, 0], 'LineWidth', 1.5)

% add parameters
set(gca, 'fontsize', 14)
ylim(yl)
xlabel('time (s)')
ylabel('power (\muV^2)')
end
function topo_plot(header, data, x_pos, x_start, map_lims)
varargin = {'maplimits' map_lims 'shading' 'interp' 'whitebk' 'on'};

% fetch data to display
x_visual = ceil((x_pos - x_start)/header.xstep);
vector = data(1, :, 1, 1, 1, x_visual);

%fetch chanlocs
chanlocs = header.chanlocs;

%parse data and chanlocs 
i=1;
for chanpos=1:size(chanlocs,2);
    vector2(i)=double(vector(chanpos));
    chanlocs2(i)=chanlocs(chanpos);
    i=i+1;
end;

topoplot(vector2,chanlocs2,varargin{:});
set(gcf,'color',[1 1 1]);
end
function [pos_x, data, sub_data] = track_peak(data, header, time_window, k, TEPs, buffer, seed)
% figure params 
figure_name = ['Peak ' TEPs.peaks{k}] ;
figure_center = TEPs.latencies(k);
span = ((1 + buffer) * TEPs.widths(k));

% set the peak timepoint manually
finish = 0;
while finish == 0
    % identify the TOI of current peak                    
    x1 = ceil((figure_center - span/2 - time_window(1)) / header.xstep);
    x2 = ceil((figure_center + span/2 - time_window(1)) / header.xstep);
    x = (figure_center - span/2) : header.xstep : (figure_center + span/2);

    % prepare data and header for visualization
    data_visual = data(1, :, :, :, :, [x1 : x2]);
    header_visual = header;
    header_visual.datasize(6) = length(data_visual);  
    header_visual.xstart = figure_center - span/2;

    % launch the figure                    
    axesHandles = [];
    fig = figure;   
    hold on          
    
    % check the length
    if length(squeeze(data_visual(1, seed, :, :, :, :))) > numel(x)
        delta = length(squeeze(data_visual(1, seed, :, :, :, :))) -  numel(x);
        data_visual = data_visual(1, :, :, :, :, 1:end-delta);
    elseif length(squeeze(data_visual(1, seed, :, :, :, :))) < numel(x)
        delta = numel(x) - length(squeeze(data_visual(1, seed, :, :, :, :)));
        data_visual(1, :, :, :, :, end + 1:end + delta) = zeros(1, delta);
    end
    
    % plot windowed data
    axesHandles = [axesHandles subplot(3, 3, [4:9])];  
    plot(x, squeeze(data_visual(1, seed, :, :, :, :)), 'LineWidth', 2)
    xlim([(figure_center - span/2), (figure_center + span/2)])
    title(figure_name, 'FontSize', 16)

    % plot the line at the center
    h = line([figure_center, figure_center], get(gca,'ylim'), 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--');   

    % plot the central topography 
    map_lims = get(axesHandles(1), 'ylim');
    axesHandles = [axesHandles subplot(3, 3, 2)];
    TEP_topoplot(header_visual, data_visual, figure_center, map_lims);

    % make the topography change with mouse movement 
    set(fig, 'WindowButtonMotionFcn', {@mouse_move, axesHandles, header_visual, data_visual});              

    % choose the peak position
    pos_x = get_position(axesHandles);             

    % update the figure
    set (fig, 'WindowButtonMotionFcn', '');
    subplot(3, 3, [4:9])
    set(h, 'XData', [pos_x, pos_x], 'YData', get(gca,'ylim'), 'LineStyle', '-');
    subplot(3, 3, 2) 
    cla(axesHandles(2))
    TEP_topoplot(header_visual, data_visual, pos_x, map_lims);

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
n_peak = ceil((pos_x - time_window(1)) / header.xstep); 

% define final length
window_length = numel((pos_x - span/2) : header.xstep : (pos_x + span/2));

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
x_out1 = ceil((pos_x - span/2 - time_window(1)) / header.xstep);
x_out2 = ceil((pos_x + span/2 - time_window(1)) / header.xstep);
sub_data = squeeze(sub_data(1, :, :, :, :, [x_out1 : x_out2]));
data = squeeze(data(1, :, :, :, :, [x_out1 : x_out2]));

function TEP_topoplot(header, data, x_pos, map_lims)
varargin = {'maplimits' map_lims 'shading' 'interp' 'whitebk' 'on'};

% fetch data to display
x_visual = ceil((x_pos - header.xstart)/header.xstep);
vector = squeeze(data(1,:,1,1,1,x_visual));

%fetch chanlocs
chanlocs = header.chanlocs;

%parse data and chanlocs 
i = 1;
for chanpos=1:size(chanlocs,2);
    vector2(i)=double(vector(chanpos));
    chanlocs2(i)=chanlocs(chanpos);
    i = i + 1;
end;

topoplot(vector2,chanlocs2,varargin{:});
set(gcf,'color',[1 1 1]);
end
function mouse_move(hObject,eventdata, axesHandles, header_visual, data_visual)

% get the position of the mouse
CP = get(hObject, 'CurrentPoint');
position = get(hObject, 'Position');
xCursor = CP(1,1)/position(1,3); % normalize
yCursor = CP(1,2)/position(1,4); % normalize

% get the position of the axes within the GUI
axesPos = get(axesHandles(1),'Position');
minx    = axesPos(1);
miny    = axesPos(2);
maxx    = minx + axesPos(3);
maxy    = miny + axesPos(4);

% check if the mouse is within the axes 
if xCursor >= minx && xCursor <= maxx && yCursor >= miny && yCursor <= maxy
    % get the cursor position within the lower axes
    currentPoint = get(axesHandles(1),'CurrentPoint');
    x = currentPoint(2,1);
    y = currentPoint(2,2);
    % update the topoplot in the uper axes
    map_lims = get(axesHandles(1), 'ylim');
    cla(axesHandles(2))
    subplot(3, 3, 2)
    TEP_topoplot(header_visual, data_visual, x, map_lims); 
    hold on
end
end
function pos_x = get_position(axesHandles)
% wait until the mouse is clicked
w = waitforbuttonpress;

% get the position of the mouse
CP = get(axesHandles(1), 'CurrentPoint');
pos_x = CP(1,1);

end
end