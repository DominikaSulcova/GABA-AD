%% RS-EEG - figures for publication
% Written by Dominika for GABA-AD project (2022)


%% params
clear all, clc

% ----- adjustable parameters -----
% dataset
prefix = 'GABA';
group = 'YC';
participant = 1:20;
medication = {'placebo' 'alprazolam'}; 
time = {'pre' 'post'};
condition = {'open' 'closed'};

% statistics
z = 1.96;
alpha = 0.05;
% --------------------------------
% navigate to the Git folder --> source of saved default files
folder_git = uigetdir(pwd, 'Choose the Git folder');

% visualization 
figure_counter = 1;

% check for colour scheme
answer = questdlg('Do you want to choose a new colour scheme?', 'Colour scheme', 'YES', 'NO', 'NO'); 
switch answer
    case 'YES'
        a = 1;
        for p = 1:6
            for c = 1:length(current)
               colours(a, :) = uisetcolor; 
               a = a + 1;
            end
        end
        save('colours.mat', 'colours'); 
    case 'NO'
    load([folder_git '\GABA_colours.mat'])
end
clear a answer

% input/output folders
folder_results = uigetdir(pwd, 'Choose the Results folder');
folder_input = [folder_results '\GABA_' group '_variables'];
folder_figures = [folder_results '\GABA_' group '_figures'];

% TOIs
TOI = struct;
TOI(1).band = 'delta'; TOI(2).band = 'theta'; TOI(3).band = 'alpha'; TOI(4).band = 'beta'; TOI(5).band = 'gamma';
TOI(1).window = [0.1, 4]; TOI(2).window = [4, 8]; TOI(3).window = [8, 13]; TOI(4).window = [13, 30]; TOI(5).window = [30, 45];
TOI(1).sign = '\delta'; TOI(2).sign = '\theta'; TOI(3).sign = '\alpha'; TOI(4).sign = '\beta'; TOI(5).sign = '\gamma';

% ROIs
ROI = struct;
ROI(1).area = 'frontal'; ROI(2).area = 'central'; ROI(3).area = 'left_temporal'; ROI(4).area = 'right_temporal'; ROI(5).area = 'occipital'; 
ROI(1).electrodes = {'Fp1' 'Fp2' 'Fz' 'F3' 'F4' 'F7' 'F8'};
ROI(2).electrodes = {'FC1' 'FC2' 'Cz' 'C1' 'C2' 'CP1' 'CP2'};
ROI(3).electrodes = {'FC5' 'T7' 'C3' 'CP5'};
ROI(4).electrodes = {'FC6' 'T8' 'C4' 'CP6'};
ROI(5).electrodes = {'P3' 'P4' 'Pz' 'P7' 'P8' 'O1' 'O2' 'Iz'};

% labels
labels = {'Fp1','Fp2','F3','F4','C3','C4','P3','P4','O1','O2',...
    'F7','F8','T7','T8','P7','P8','Fz','Cz','Pz','Iz','FC1','FC2',...
    'CP1','CP2','FC5','FC6','CP5','CP6','P5', 'P6', 'C1','C2'};

%% higher plot spectra by medication 
% ----- section input -----
region = 'central';
crop_window = [1, 30]; 
% -------------------------
% load data and header
load([folder_input '\rsEEG_data_high.mat']);
data = data_high; clear data_high
load([folder_git '\header_high.mat']);
header = header_high; clear header_high

% determinde region
for r = 1:length(ROI)
    if strcmp(ROI(r).area, region)
        region_n = r;
    end
end

% determine x parameters
x_start = ceil(crop_window(1) - header.xstart)/header.xstep;
x_end = ceil(crop_window(2) - header.xstart)/header.xstep;
x = [crop_window(1) : header.xstep : crop_window(2)];

%%
% loop through medications
for m = 1:length(medication)
    % subset the data
    data_bl = squeeze(mean(data(m, 1, :, :, region_n, x_start:x_end), 4));
    SEM_bl =  squeeze(std(data(m, 1, :, :, region_n, x_start:x_end), 0, 4)/sqrt(size(data, 4)));
    data_post = squeeze(mean(data(m, 2, :, :, region_n, x_start:x_end), 4));
    SEM_post =  squeeze(std(data(m, 2, :, :, region_n, x_start:x_end), 0, 4)/sqrt(size(data, 4)));
    
    % launch the figure
    fig = figure(figure_counter);
    hold on
    
    % set up limits
    ylim([-1.6, -0.2])
    xlim([-1, 32])
    
    % plot lines between frequency bands (= TOIs)
    for h = 1:length(TOI)-2
        line([TOI(h).window(2), TOI(h).window(2)], [-1.6, -0.2], 'Color', [0.9 0.9 0.9], 'LineWidth', 0.8);
        hold on
    end
    
    % add names of frequency bands
    for h = 1:length(TOI) - 1
        text((TOI(h).window(1) + TOI(h).window(2))/2, -1.45, TOI(h).sign, 'fontsize', 16, 'color', [0.6 0.6 0.6])
        hold on
    end
    
    % plot data
    plot_data(x, data_bl, SEM_bl, [0 0 0])
    plot_data(x, data_post, SEM_post, colours((m-1)*2 + 2, :))   
    
    % save figure
    figure_name = sprintf('rsEEG_spectra_%s', medication{m});
    savefig([folder_figures '\' figure_name '.fig'])
    saveas(fig, [folder_figures '\' figure_name '.svg'], 'svg')
    
    % updatefigure counter
    figure_counter = figure_counter + 1;
       
end
%%

clear r m x region_n data header data_bl SEM_bl data_post SEM_post fig
clear linestyles shades yl h c F P
clear region freq crop_window

%% functions
function plot_data(x, data, SEM, colour)
    % params 
    linestyles = {'-', ':'};
    shades = [0.3 0.15];
    
    % plot SEM
    for c = 1:size(data, 1)
        F(c) = fill([x fliplr(x)],[log10(data(c, :) + SEM(c, :)) log10(fliplr(data(c, :) - SEM(c, :)))], ...
        colour, 'FaceAlpha', shades(c), 'linestyle', 'none');
        hold on
    end
    
    % plot data
    for c = 1:size(data, 1)
        P(c) = plot(x, log10(data(c, :)), 'Color', colour, 'LineWidth', 2.5, 'LineStyle', linestyles{c});
        hold on
    end
end
