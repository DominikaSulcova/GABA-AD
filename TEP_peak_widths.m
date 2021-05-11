%% script that extracts default peak widths and latencies
% 1) loads grand average baseline dataaset
%% parameters
% dataset
dataset = 'avg merged_subj all pre';
electrode = {'Cz' 'target'};
load('labels.mat')

% TOI
time_window = [-0.005, 0.4];

%% prepare the data
load([dataset '.mat'])
load([dataset '.lw6'], '-mat')

% index Cz and target
Cz = find(contains(labels, 'Cz')); 
target = find(contains(labels, 'target')); 

% extract data from channels of interest
Cz_data = double(squeeze(data(1, Cz, :, :, :, [((time_window(1)-header.xstart)/header.xstep):((time_window(2)-header.xstart)/header.xstep)+1])));
target_data = double(squeeze(data(1, target, :, :, :, [((time_window(1)-header.xstart)/header.xstep):((time_window(2)-header.xstart)/header.xstep)+1])));

x = [time_window(1):header.xstep:time_window(2)];

% check data by plotting
figure(1)
plot(x, Cz_data, 'b')
xlim([time_window(1), time_window(2)])
hold on
plot(x, target_data, 'r')
%% find peak widths
% prepare the final table
peak_widths = table;
peak_widths.electrode = {'Cz'; 'Cz'; 'Cz'; 'Cz'; 'Cz'; 'Cz'; ...
    'target'; 'target'; 'target'; 'target'; 'target'; 'target'};

% Cz electrode
[peaks_pos, locs_pos, width_pos, promi_pos] = findpeaks(Cz_data, 'Annotate','extents','MinPeakDistance', 0.02, 'MinPeakProminence', 0.5);
[peaks_neg, locs_neg, width_neg, promi_neg] = findpeaks(-Cz_data, 'Annotate','extents','MinPeakDistance', 0.02, 'MinPeakProminence', 0.4);

locs = ([locs_pos', locs_neg'] * header.xstep) + time_window(1);
locs_Cz = locs([4, 1, 5, 2, 6, 3]);

widths = [width_pos', width_neg'] * header.xstep;
widths_Cz = widths([4, 1, 5, 2, 6, 3]);

% target electrodes
[peaks_pos, locs_pos, width_pos, promi_pos] = findpeaks(target_data, 'Annotate','extents','MinPeakDistance', 0.02, 'MinPeakProminence', 0.5);
[peaks_neg, locs_neg, width_neg, promi_neg] = findpeaks(-target_data, 'Annotate','extents','MinPeakDistance', 0.02, 'MinPeakProminence', 0.4);

locs = ([locs_pos', locs_neg'] * header.xstep) + time_window(1);
locs_target = locs([4, 1, 5, 2, 6, 3]);

widths = [width_pos', width_neg'] * header.xstep;
widths_target = widths([4, 1, 5, 2, 6, 3]);

% deal with the peak widths within the positive complex P30-N45-P60
width_complex = widths_target(4);

y = (peaks_pos(2) - promi_pos(2)/2);
h = ones(1, length(x)) * y;                                                 % h is a x-long vector of monotonous y-value corresponding to the height of the half of the prominence of the positive complex 
[x_h,y_h] = intersections(x,target_data,x,h,1);                             % calculate intersections of h and target TEP curve
h_N45 = ones(1, length(x)) * (-peaks_neg(2) + promi_neg(2)/2);              % x-long vector of monotonous y-value corresponding to the height of the half of the prominence of N45
[x_N45,y_N45] = intersections(x,target_data,x,h_N45,1);                     % calculate intersections of h_N45 and target TEP curve

widths_target(2) = x_N45(2) - x_h(1);
widths_target(4) = x_h(2) - x_N45(3);

floor(widths_target(2) + widths_target(3) + widths_target(4)) == floor(width_complex)

% plot the division of P30-N45-P60
figure(3)
plot(x, target_data, 'b', 'LineWidth', 2)
hold on
xlim([-0.005, 0.15])
plot(x, h, '--k')
% plot borders of the interval of positive complex P30-N45-P60
line([x_h(1), x_h(1)], get(gca,'ylim'), 'LineStyle', '--', 'Color', [0, 0, 0])
line([x_h(2), x_h(2)], get(gca,'ylim'), 'LineStyle', '--', 'Color', [0, 0, 0])
% plot borders of N45 interval
line([x_N45(2), x_N45(2)], get(gca,'ylim'), 'LineStyle', '--', 'Color', [1 0 0])
line([x_N45(3), x_N45(3)], get(gca,'ylim'), 'LineStyle', '--', 'Color', [1 0 0])
% highlight widths of peaks P30, N45, and P60
line([x_h(1), x_N45(2)], [y, y], 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2)
line([x_N45(2), x_N45(3)], [y, y], 'Color', [1 0 0], 'LineWidth', 2)
line([x_N45(3), x_h(2)], [y, y], 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2)
% highlight crossections
plot(x_h([1, 2]), y_h([1, 2]), 'ro')
plot(x_N45([2, 3]), y_N45([2, 3]), 'ro')
plot(x_N45([2, 3]), y_h([1, 2]), 'ro')
hold off

% fill in the table
locs = [locs_Cz, locs_target];
widths = ceil([widths_Cz, widths_target] * 1000);
peak_widths.latency = locs';
peak_widths.width = widths';

save('peak_widths.mat', 'peak_widths')

% plot the peak widths
figure(3)
findpeaks(Cz_data, x, 'Annotate','extents','MinPeakDistance', 0.02, 'MinPeakProminence', 0.5);
hold on
plot(x, Cz_data, 'r', 'LineWidth', 2)
findpeaks(-Cz_data, x,'Annotate','extents','MinPeakDistance', 0.02, 'MinPeakProminence', 0.4)
plot(x, -Cz_data, '--r', 'LineWidth', 2)
title('Peak widths - Cz electrode')
hold off

figure(4)
findpeaks(target_data, x, 'Annotate','extents','MinPeakDistance', 0.02, 'MinPeakProminence', 0.5);
hold on
plot(x, target_data, 'r', 'LineWidth', 2)
findpeaks(-target_data, x,'Annotate','extents','MinPeakDistance', 0.02, 'MinPeakProminence', 0.4)
plot(x, -target_data, '--r', 'LineWidth', 2)
title('Peak widths - target electrodes')
hold off



