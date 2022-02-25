%% GABA-AD: TEPs - TOPOGRAPHICAL ANALYSIS
% Written by Dominika for GABA-AD project (2022)
% 
% Colection of scripts to visualize the outcome of the topographical analysis performed in Ragu 
%   --> figures are saved in a folder 'GABA_YC_figures'
% 
% 1) load the data
% 

%% parameters
clear all; clc

% ----- adjustable parameters -----
% dataset
prefix = 'GABA';
group = 'YC';
participant = 1:20;
medication = {'placebo' 'alprazolam'}; 
time = {'pre' 'post'};
stimulus = {'M1 - CS' 'M1 - TS' 'AG' 'M1 - ppTMS'};
peaks_M1 = {'N17' 'P30' 'N45' 'P60' 'N100' 'P180'};
peaks_AG = {'P25' 'N40' 'P50' 'N75' 'N100' 'P180'};

% visualization (in ms)
time_window = [-50, 300];
analysis_window = [10, 300];
x_delta = 0.5;
shade = 0.2;

% statistics
z = 1.96;
alpha = 0.05;
% --------------------------------

% navigate to the Git folder --> source of saved default files
folder_git = uigetdir(pwd, 'Choose the Git folder');

% visualization 
figure_counter = 1;
x = time_window(1):x_delta:time_window(2);

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
folder_input = [folder_results '\GABA_' group '_microstates'];
folder_figures = [folder_results '\GABA_' group '_figures'];

%% 1) LOAD DATA
data = struct;

% baseline - target
load([folder_input '\GABA_' group '_baseline_target.mat'])
data.target = rd;

% baseline - intensity
load([folder_input '\GABA_' group '_baseline_intensity.mat'])
data.intensity = rd;

% medication - M1 TS
load([folder_input '\GABA_' group '_medication_M1-TS.mat'])
data.M1_TS = rd;

% medication - M1 CS
load([folder_input '\GABA_' group '_medication_M1-CS.mat'])
data.M1_CS = rd;

% medication - AG
load([folder_input '\GABA_' group '_medication_AG.mat'])
data.AG = rd;

% set up condiitions
condition = fieldnames(data);
clear rd

%% ) TANOVA OUTCOME 
for d = 1:numel(condition)
    % extract comparison names
    statement = ['comps = {data.' condition{d} '.strF1, data.' condition{d} '.strF2, [data.' condition{d} '.strF1 '' * '' data.' condition{d} '.strF2]};'];
    eval(statement)
    
    % TANOVA p value
    for a = 1:length(comps)
        % load data
        statement = ['data_visual(a, :) = squeeze(data.' condition{d} '.PTanova(1, a+1, :, 1));'];
        eval(statement)
        
        % define intervals of significance
        signif_05 = data_visual < 0.05;
        signif_01 = data_visual < 0.01;
        
        % plot the p-value timecourse
        % launch the figure
        fig = figure(figure_counter);
        hold on
        
        % shade intervals of significance
        I(1) = area(x, signif_05);
        I(1).FaceColor = [1 0.73 0.73];
        I(1).EdgeColor = 'none';
        I(2) = area(x, signif_01);
        I(2).FaceColor = [1 0.44 0.44];
        I(2).EdgeColor = 'none';
        
        % plot p
        P = plot(x, data_visual, 'Color', [0 0 0], 'LineWidth', 3);
        
        % other parameters
        title([condition{d} ' TANOVA: factor ''' comps{a} ''''])
        xlabel('time (ms)')
        ylabel('probability')
        set(gca, 'FontSize', 18)
        xlim([analysis_window(1) + 1, analysis_window(2)] )    
        hold off

        % change figure size
        fig.Position = [500 500 750 300];
        
    end   

end
clear d a statement signif_05 signif_01 fig I 


