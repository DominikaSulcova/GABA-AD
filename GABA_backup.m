%% DISS
% m = 1; t = 1; s = 1, p = 1; i = 1;
% normalize data by GMFP
for m = 1:length(medication)
    for t = 1:length(time)
        for s = 1:length(stimulus) 
            for p = 1:length(participant)
                for i = 1:size(GABA_data, 6)
                    % devide data at each time point by GMFP
                    GABA_data_norm(m, t, s, p, :, i) = squeeze(GABA_data(m, t, s, p, 1:30, i)) / GABA_GMFP_subject(m, t, s, p, i);
                end
            end
        end
    end
end
clear m t s p i

% ----- baseline -----
% select baseline data of both sessions
for s = 1:length(stimulus)
    counter = 1;
    for p = 1:length(participant)
        for m = 1:length(medication)
            data_temp(s, counter, :, :) = squeeze(GABA_data_norm(m, 1, s, p, :, :));
            counter = counter + 1;
        end
    end
end
clear s p m counter

% average across subjects
for s = 1:length(stimulus)
    data(s, :, :) = squeeze(mean(data_temp(s, :, :, :), 2));
end
clear s data_temp

% calculate DISS and spatial correlation C
GABA_DISS = struct;
for i = 1:size(GABA_data_norm, 6)
    % between TS and CS stimuli
    diff = squeeze(data(2, :, i) - data(1, :, i));
    GABA_DISS.baseline.stimulus(1, i) = sqrt(mean(diff.^2));
    GABA_DISS.baseline.stimulus(2, i) = 1 - (GABA_DISS.baseline.stimulus(1, i)^2)/2;

    % between ppTMS and TS stimuli
    diff = squeeze(data(3, :, i) - data(2, :, i));
    GABA_DISS.baseline.SICI(1, i) = sqrt(mean(diff.^2));
    GABA_DISS.baseline.SICI(2, i) = 1 - (GABA_DISS.baseline.SICI(1, i)^2)/2;
end
clear p i diff

% plot baseline TS-CS
data_visual = GABA_DISS.baseline.stimulus;
fig = figure(figure_counter);
hold on
title('Global dissimilarity: baseline TS - CS', 'FontSize', 16, 'FontWeight', 'bold')
set(gca, 'fontsize', 12); xlabel('time (s)')
ax = gca; ax.XColor = [0.5 0.5 0.5];
plot(x, data_visual(1, :), 'Color', colours(2, :), 'LineWidth', 2.5)
fill([x fliplr(x)],[data_visual(2, :) zeros(1, length(data_visual(2, :)))], ...
    colours(1, :) , 'facealpha', 0.5, 'linestyle', 'none');
line(time_window, [0 0], 'Color', [0, 0, 0], 'LineWidth', 0.5)
line(time_window, [1 1], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5)
line([0, 0], [-0.75 1.75], 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5) 
ylim([-0.75 1.75]); 
lgd = legend({'DISS index' 'spatial correlation'}, 'Location', 'southeast');
lgd.FontSize = 12;
hold off
figure_name = ['DISS_' target '_bl_stimulus'];
savefig([folderpath '\' filename '_figures\' figure_name '.fig'])
saveas(fig, [folderpath '\' filename '_figures\' figure_name '.png'])
figure_counter = figure_counter + 1;
clear data_visual fig ax lgd figure_name

% plot baseline ppTMS-TS (SICI)
data_visual = GABA_DISS.baseline.SICI;
fig = figure(figure_counter);
hold on
title('Global dissimilarity: baseline SICI', 'FontSize', 16, 'FontWeight', 'bold')
set(gca, 'fontsize', 12); xlabel('time (s)')
ax = gca; ax.XColor = [0.5 0.5 0.5];
plot(x, data_visual(1, :), 'Color', colours(6, :), 'LineWidth', 2.5)
fill([x fliplr(x)],[data_visual(2, :) zeros(1, length(data_visual(2, :)))], ...
    colours(5, :) , 'facealpha', 0.5, 'linestyle', 'none');
line(time_window, [0 0], 'Color', [0, 0, 0], 'LineWidth', 0.5)
line(time_window, [1 1], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5)
line([0, 0], [-0.75 1.75], 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5) 
ylim([-0.75 1.75]); 
lgd = legend({'DISS index' 'spatial correlation'}, 'Location', 'southeast');
lgd.FontSize = 12;
hold off
figure_name = ['DISS_' target '_bl_SICI'];
savefig([folderpath '\' filename '_figures\' figure_name '.fig'])
saveas(fig, [folderpath '\' filename '_figures\' figure_name '.png'])
figure_counter = figure_counter + 1;
clear data_visual fig ax lgd figure_name

% ----- pre-post medication -----
% calculate DISS and spatial correlation C
for m = 1:length(medication)
    for s = 1:length(stimulus)
        for i = 1:size(GABA_data_norm, 6)
            diff = mean(squeeze(GABA_data_norm(m, 2, s, :, :, i)), 1) - mean(squeeze(GABA_data_norm(m, 1, s, :, :, i)), 1);
            GABA_DISS.medication(1, m, s, i) = sqrt(mean(diff.^2));
            GABA_DISS.medication(2, m, s, i) = 1 - (GABA_DISS.medication(1, m, s, i)^2)/2;
        end
    end
end
clear m s i diff

% plot 
for m = 1:length(medication)
    for s = 1:length(stimulus) 
        data_visual = squeeze(GABA_DISS.medication(:, m, s, :));
        fig = figure(figure_counter);
        hold on
        title(sprintf('Global dissimilarity: %s, %s', medication{m}, stimulus{s}), 'FontSize', 16, 'FontWeight', 'bold')
        set(gca, 'fontsize', 12); xlabel('time (s)')
        ax = gca; ax.XColor = [0.5 0.5 0.5];
        plot(x, data_visual(1, :), 'Color', colours((m-1)*2 + 2, :), 'LineWidth', 2.5)
        fill([x fliplr(x)],[data_visual(2, :) zeros(1, length(data_visual(2, :)))], ...
            colours((m-1)*2 + 1, :) , 'facealpha', 0.5, 'linestyle', 'none');
        line(time_window, [0 0], 'Color', [0, 0, 0], 'LineWidth', 0.5)
        line(time_window, [1 1], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5)
        line([0, 0], [-0.75 1.75], 'LineStyle', '--', 'Color', [0, 0, 0], 'LineWidth', 2.5) 
        ylim([-0.75 1.75]); 
        lgd = legend({'DISS index' 'spatial correlation'}, 'Location', 'southeast');
        lgd.FontSize = 12;
        hold off
        figure_name = ['DISS_' target '_' medication{m} '_' stimulus{s}];
        savefig([folderpath '\' filename '_figures\' figure_name '.fig'])
        saveas(fig, [folderpath '\' filename '_figures\' figure_name '.png'])
        figure_counter = figure_counter + 1;
    end
end
clear data_visual fig ax lgd figure_name

%% peak widths
% calculate narrow window parameters
x_start_narrow = (0.01 - time_window(1))/xstep;
x_end_narrow = (0.25 - time_window(1))/xstep;
x_narrow = [0.01:xstep:0.25];

% loop through stimuli
for s = 1:length(stimulus)
    % average data across conditions
    for i = 1:length(x_narrow)
        data_visual(i) = mean(GABA_GMFP(:, :, s, x_start_narrow + (i - 1)), 'all');
    end
        
    % identify peak widths
    [P, L, GABA_peaks(s).widths, R] = findpeaks(data_visual, 'Annotate','extents', 'WidthReference', 'halfheight');
    GABA_peaks(s).widths = ceil(GABA_peaks(s).widths)* header.xstep;
    
    % plot figure
    fig = figure(figure_counter);
    hold on
    findpeaks(data_visual, x_narrow, 'Annotate', 'extents', 'WidthReference', 'halfheight');
    set(gca, 'fontsize', 14)
    xlabel('time(s)')
    ylabel('power (\muV^2)')
    grid off
    
    % add width denotation
    for k = 1:length(GABA_peaks(s).widths)
        if k == length(GABA_peaks(s).widths)
            text(GABA_peaks(s).latencies(k) - 0.005, -0.25, sprintf('%1.0f ms', GABA_peaks(s).widths(k)*1000), 'Color', [0.93 0.69 0.13], 'FontSize', 14)
        else
            text(GABA_peaks(s).latencies(k) - 0.005, -0.25, sprintf('%1.0f', GABA_peaks(s).widths(k)*1000), 'Color', [0.93 0.69 0.13], 'FontSize', 14)
        end
    end
    xlim([-0.005, 0.26])
    
    % save figure, update
    savefig(['GABA_' target '_widths'])
    saveas(fig, ['GABA_' target '_widths.png'])
    figure_counter = figure_counter + 1;
    end
clear s i

%%% ONLY FOR AG %%%
GABA_peaks(s).widths = [GABA_peaks(s).widths(1:2), GABA_peaks(s).widths(2), GABA_peaks(s).widths(3:5)];
%%%%%%%%%%%%%%%%%%%

% append new variables to the general MATLAB file
save(filename, 'GABA_TEP', '-append');
clear x_start_narrow x_end_narrow x_narrow s k data fig P L R 

%% functions
function [peak_x, peak_y] = gfp_peaks(y, time_window, xstep, varargin)
% check whether to plot labels (default)
if ~isempty(varargin)
    a = find(strcmpi(varargin, 'max_peaks'));
    if ~isempty(a)
        max_peaks = varargin{a + 1};
    end
end

% find peaks 
[pks, locs] = findpeaks(y, 'MinPeakDistance', 5, 'MinPeakProminence', 0.01);
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
function [p1, p2, f] = gfp_plot_diff(x, y, time_window, colours)
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