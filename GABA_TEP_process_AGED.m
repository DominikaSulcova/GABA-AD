%% GABA-AD: TMS-EVOKED POTENTIALS
% Written by Dominika for GABA-AD project (2021)

%% parameters
clear all; clc

% dataset
group = {'MCI' 'MCI-CTRL'};
participant = 1:20;
medication = {'placebo' 'alprazolam'}; 
time = {'pre' 'post'};
prefix = 'avg bl icfilt ica26 crop but fft-notchfilt prefilt prea28 visual reref ds art-sup bl ep dc chan-select';

% target & stimulus
target = questdlg('Choose target:', 'Stimulated cortex', 'M1', 'AG', 'M1');
switch target
    case 'M1'
        stimulus = {'CS' 'TS' 'ppTMS'};
    case 'AG'
        stimulus = {''};
end

% electrode labels
load('labels.mat')
switch target
    case 'M1'
        labels.CS = labels_CS; labels.TS = labels_TS;
    case 'AG'
        labels = labels_AG;        
end
clear labels_CS labels_TS labels_AG

% visualization
load([prefix ' ' group{1} num2str(participant(1)) ' ' medication{1} ' ' target ' ' time{1} '.lw6'], '-mat')
xstep = header.xstep; 
time_window = [-0.05, 0.3];
x = [time_window(1):xstep:time_window(2)];
x_start = (time_window(1) - header.xstart)/xstep;
x_end = (time_window(2) - header.xstart)/xstep;
load('colours.mat'); colours = colours([2 1 4 6],:);
alpha = 0.2;
figure_counter = 1;
clear header

% stats
z = 1.96;




