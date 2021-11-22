%% GABA-AD: TMS-EVOKED POTENTIALS 
% Written by Dominika for GABA-AD project (2020-21)
% 
% Colection of scripts to visualy explore extracted TMS-EEG variables: 
%   --> output datasets are saved in a folder 'GABA_YC_statistics'
%   --> figures are saved in a folder 'GABA_YC_figures'
% 
% 1) load the data


%% parameters
clear all; clc

% ----- adjustable parameters -----
% dataset
prefix = 'GABA';
group = 'YC';
participant = 1:20;
medication = {'placebo' 'alprazolam'}; 
time = {'pre' 'post'};
stimulus = {'CS' 'TS' 'ppTMS'};
input_file_M1 = [prefix '_' group '_M1_TEPs.mat'];
input_file_AG = [prefix '_' group '_AG_TEPs.mat'];

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

% create output folders
folder_results = uigetdir(pwd, 'Choose the Results folder');
folder_input = [folder_results '\GABA_' group '_variables'];
folder_output = [folder_results '\GABA_' group '_statistics'];
folder_figures = [folder_results '\GABA_' group '_figures'];
if ~exist(folder_output)
    mkdir(folder_output)
end
if ~exist(folder_figures)
    mkdir(folder_figures)  
end
clear folder_results

%% 1) EFFECT OF ALPRAZOLAM - TEPs: load data
% parameters
load([folder_output '\GABA_YC_results.mat'])

% TEPs from M1
load([folder_input '\' input_file_M1], 'GABA_TEP_peaks', 'GABA_TEP_default');
data_M1 = squeeze(GABA_TEP_peaks.change(:, :, 1, 2:6, [2 3]));
peaks_M1 = GABA_TEP_default.peak(2:6); 
sprintf('M1 - datasize: %s', num2str(size(data_M1)))

% TEPs from AG
load([folder_input '\' input_file_AG], 'GABA_TEP_peaks', 'GABA_TEP_default');
data_AG = squeeze(GABA_TEP_peaks.change(:, :, :, [2 3]));
peaks_AG = GABA_TEP_default.peak; 
sprintf('AG - datasize: %s', num2str(size(data_AG)))

clear GABA_TEP_peaks GABA_TEP_default

%% 2) EFFECT OF ALPRAZOLAM - TEPs: distribution
% M1: hostogram - peak amplitude
fig = figure(figure_counter);
for p = 1:length(peaks_M1)
    for m = 1:length(medication)        
        % launch the plot
        subplot(length(medication), length(peaks_M1), (m-1)*length(peaks_M1) + p)
        hold on
        
        % plot the histogram with a normal distribution fit
        h = histfit(data_M1(:, m, p, 1), 10, 'kernel');
        h(1).FaceColor = colours((m-1)*2 + 1, :); h(1).EdgeColor = [1 1 1];
        h(2).Color = [0 0 0]; h(2).LineWidth = 1;
        
        % get limits
        yl = get(gca, 'ylim'); xl = get(gca, 'xlim');
        
        % add mean of the distribution
        mu = mean(data_M1(:, m, p, 1));
        text(xl(1), yl(2)*-0.1, sprintf('mu = %1.3f', mu), 'FontSize', 10, 'FontWeight', 'bold')
        
        % add skewness and kurtosis
        s = skewness(data_M1(:, m, p, 1)); k = kurtosis(data_M1(:, m, p, 1));
        text(xl(1), yl(2)*-0.2, sprintf('%1.3f, %1.3f', s, k), 'FontSize', 10)
            
        % set annotation
        if p == 1
            text(xl(1) - (xl(2) - xl(1)), yl(2)*0.7, medication{m}, 'FontSize', 10)
            clear yl xl
        end
                       
        % set title
        if m == 1
            title(peaks_M1{p})
        end
    end
end   
clear p m h yl xl mu s k

% parameters
fig.Position = [500 400 750 500];
sgtitle('M1 - change in peak amplitude')

% update counter
figure_counter = figure_counter + 1;

% M1: QQ plot - peak amplitude
fig = figure(figure_counter);
for p = 1:length(peaks_M1)
    for m = 1:length(medication)        
        % launch the plot
        subplot(length(medication), length(peaks_M1), (m-1)*length(peaks_M1) + p)
        hold on
        
        % plot the histogram with a normal distribution fit
        q = qqplot(data_M1(:, m, p, 1));
        title(''); xlabel(''); ylabel('')
        
        % get limits
        yl = get(gca, 'ylim'); xl = get(gca, 'xlim');
      
        % set annotation
        if p == 1
            text(xl(1) - (xl(2) - xl(1))*0.7, yl(2)*0.7, medication{m}, 'FontSize', 10)
            clear yl xl
        end
                       
        % set title
        if m == 1
            title(peaks_M1{p})
        end
    end
end   
clear p m q yl xl

% parameters
fig.Position = [200 400 1000 500];
sgtitle('M1 - change in peak amplitude')

% update counter
figure_counter = figure_counter + 1;

% AG: hostogram - peak amplitude
fig = figure(figure_counter);
for p = 1:length(peaks_AG)
    for m = 1:length(medication)        
        % launch the plot
        subplot(length(medication), length(peaks_AG), (m-1)*length(peaks_AG) + p)
        hold on
        
        % plot the histogram with a normal distribution fit
        h = histfit(data_AG(:, m, p, 1), 10, 'kernel');
        h(1).FaceColor = colours((m-1)*2 + 1, :); h(1).EdgeColor = [1 1 1];
        h(2).Color = [0 0 0]; h(2).LineWidth = 1;
        
        % get limits
        yl = get(gca, 'ylim'); xl = get(gca, 'xlim');
        
        % add mean of the distribution
        mu = mean(data_AG(:, m, p, 1));
        text(xl(1), yl(2)*-0.1, sprintf('mu = %1.3f', mu), 'FontSize', 10, 'FontWeight', 'bold')
        
        % add skewness and kurtosis
        s = skewness(data_AG(:, m, p, 1)); k = kurtosis(data_AG(:, m, p, 1));
        text(xl(1), yl(2)*-0.2, sprintf('%1.3f, %1.3f', s, k), 'FontSize', 10)
            
        % set annotation
        if p == 1
            text(xl(1) - (xl(2) - xl(1)), yl(2)*0.7, medication{m}, 'FontSize', 10)
            clear yl xl
        end
                       
        % set title
        if m == 1
            title(peaks_AG{p})
        end
    end
end   
clear p m h yl xl mu s k

% parameters
fig.Position = [500 400 750 500];
sgtitle('AG - change in peak amplitude')

% update counter
figure_counter = figure_counter + 1;

% AG: QQ plot - peak amplitude
fig = figure(figure_counter);
for p = 1:length(peaks_AG)
    for m = 1:length(medication)        
        % launch the plot
        subplot(length(medication), length(peaks_AG), (m-1)*length(peaks_AG) + p)
        hold on
        
        % plot the histogram with a normal distribution fit
        q = qqplot(data_AG(:, m, p, 1));
        title(''); xlabel(''); ylabel('')
        
        % get limits
        yl = get(gca, 'ylim'); xl = get(gca, 'xlim');
      
        % set annotation
        if p == 1
            text(xl(1) - (xl(2) - xl(1))*0.7, yl(2)*0.7, medication{m}, 'FontSize', 10)
            clear yl xl
        end
                       
        % set title
        if m == 1
            title(peaks_AG{p})
        end
    end
end   
clear p m q yl xl 

% parameters
fig.Position = [200 400 1000 500];
sgtitle('AG - change in peak amplitude')

% update counter
figure_counter = figure_counter + 1;
clear fig

%% 3) EFFECT OF ALPRAZOLAM - TEPs: square root
% add a constant to make all entries positive
c = [floor(min(data_M1, [], 'all')) floor(min(data_AG, [], 'all'))];
data_M1_sqrt = data_M1 - c(1);
data_AG_sqrt = data_AG - c(2);
clear c

% % reflect negatively skewed data
% for p = find(skewness == -1)
%     for m = 1:length(medication)
%         for t = 1:length(time)
%             for s = 1:length(stimulus)
%                 GABA_TEP_peaks.amplitude_peak_log(m, t, s, :, p) = max(GABA_TEP_peaks.amplitude_peak_log(m, t, s, :, p)) - GABA_TEP_peaks.amplitude_peak_log(m, t, s, :, p) + 1;
%             end
%         end
%     end
% end
% clear p m t s 

% slog transform the data
data_M1_sqrt = sqrt(data_M1_sqrt);
data_AG_sqrt = sqrt(data_AG_sqrt);

% replot
% M1: histogram - peak amplitude
fig = figure(figure_counter);
for p = 1:length(peaks_M1)
    for m = 1:length(medication)        
        % launch the plot
        subplot(length(medication), length(peaks_M1), (m-1)*length(peaks_M1) + p)
        hold on
        
        % plot the histogram with a normal distribution fit
        h = histfit(data_M1_sqrt(:, m, p, 1), 10, 'kernel');
        h(1).FaceColor = colours((m-1)*2 + 1, :); h(1).EdgeColor = [1 1 1];
        h(2).Color = [0 0 0]; h(2).LineWidth = 1;
        
        % get limits
        yl = get(gca, 'ylim'); xl = get(gca, 'xlim');
        
        % add mean of the distribution
        mu = mean(data_M1_sqrt(:, m, p, 1));
        text(xl(1), yl(2)*-0.1, sprintf('mu = %1.3f', mu), 'FontSize', 10, 'FontWeight', 'bold')
        
        % add skewness and kurtosis
        s = skewness(data_M1_sqrt(:, m, p, 1)); k = kurtosis(data_M1_sqrt(:, m, p, 1));
        text(xl(1), yl(2)*-0.2, sprintf('%1.3f, %1.3f', s, k), 'FontSize', 10)
            
        % set annotation
        if p == 1
            text(xl(1) - (xl(2) - xl(1)), yl(2)*0.7, medication{m}, 'FontSize', 10)
            clear yl xl
        end
                       
        % set title
        if m == 1
            title(peaks_M1{p})
        end
    end
end   
clear p m h yl xl mu s k

% parameters
fig.Position = [500 400 750 500];
sgtitle('M1 - change in peak amplitude (square root)')

% update counter
figure_counter = figure_counter + 1;

% AG: histogram - peak amplitude
fig = figure(figure_counter);
for p = 1:length(peaks_AG)
    for m = 1:length(medication)        
        % launch the plot
        subplot(length(medication), length(peaks_M1), (m-1)*length(peaks_M1) + p)
        hold on
        
        % plot the histogram with a normal distribution fit
        h = histfit(data_AG_sqrt(:, m, p, 1), 10, 'kernel');
        h(1).FaceColor = colours((m-1)*2 + 1, :); h(1).EdgeColor = [1 1 1];
        h(2).Color = [0 0 0]; h(2).LineWidth = 1;
        
        % get limits
        yl = get(gca, 'ylim'); xl = get(gca, 'xlim');
        
        % add mean of the distribution
        mu = mean(data_AG_sqrt(:, m, p, 1));
        text(xl(1), yl(2)*-0.1, sprintf('mu = %1.3f', mu), 'FontSize', 10, 'FontWeight', 'bold')
        
        % add skewness and kurtosis
        s = skewness(data_AG_sqrt(:, m, p, 1)); k = kurtosis(data_AG_sqrt(:, m, p, 1));
        text(xl(1), yl(2)*-0.2, sprintf('%1.3f, %1.3f', s, k), 'FontSize', 10)
            
        % set annotation
        if p == 1
            text(xl(1) - (xl(2) - xl(1)), yl(2)*0.7, medication{m}, 'FontSize', 10)
            clear yl xl
        end
                       
        % set title
        if m == 1
            title(peaks_AG{p})
        end
    end
end   
clear p m h yl xl mu s k

% parameters
fig.Position = [500 400 750 500];
sgtitle('AG - change in peak amplitude (square root)')

% update counter
figure_counter = figure_counter + 1;
clear fig

%% 4) EFFECT OF ALPRAZOLAM - TEPs: export for R
GABA_YC_medication_TEP = table;
row_counter = 1;
for p = 1:length(participant) 
    for m = 1:length(medication)  
        % independent variables
        GABA_YC_medication_TEP.subject(row_counter) = participant(p);
        GABA_YC_medication_TEP.medication(row_counter) = medication(m); 
                
        % covariate
        GABA_YC_medication_TEP.rMT(row_counter) = GABA_YC_results.rmt(m).change(p);  
                
        % dependent variables - M1 
        for k = 1:length(peaks_M1)             
            statement = ['GABA_YC_medication_TEP.M1_amp_' peaks_M1{k} '(row_counter) = data_M1(p, m, k, 1);'];
            eval(statement)
            statement = ['GABA_YC_medication_TEP.M1_lat_' peaks_M1{k} '(row_counter) = data_M1(p, m, k, 2);'];
            eval(statement)                
        end
        
        % dependent variables - AG 
        for k = 1:length(peaks_AG)             
            statement = ['GABA_YC_medication_TEP.AG_amp_' peaks_AG{k} '(row_counter) = data_AG(p, m, k, 1);'];
            eval(statement)
            statement = ['GABA_YC_medication_TEP.AG_lat_' peaks_AG{k} '(row_counter) = data_AG(p, m, k, 2);'];
            eval(statement)                
        end
                    
        % update the counter
        row_counter = row_counter + 1;
    end
end
clear m t s p k row_counter
writetable(GABA_YC_medication_TEP, [folder_output '\GABA_YC_medication_TEP.csv'])

%% 5) EFFECT OF ALPRAZOLAM - RS-EEG: load the data
for m = 1:length(medication)
    data_rsEEG(:, m, 1) = GABA_YC_results.rsEEG(m).sigma.change';  
    data_rsEEG(:, m, 2) = GABA_YC_results.rsEEG(m).delta.change'; 
    data_rsEEG(:, m, 3) = GABA_YC_results.rsEEG(m).AAC.change'; 
    data_rsEEG(:, m, 4) = GABA_YC_results.rsEEG(m).SE.change.closed(:, 1)'; 
end
clear m
sprintf('M1 - datasize: %s', num2str(size(data_rsEEG)))
DV_rsEEG = {'sigma' 'delta' 'AAC' 'SE'};

%% 6) EFFECT OF ALPRAZOLAM - RS-EEG: distribution
% histogram 
fig = figure(figure_counter);
for p = 1:length(DV_rsEEG)
    for m = 1:length(medication)        
        % launch the plot
        subplot(length(medication), length(DV_rsEEG), (m-1)*length(DV_rsEEG) + p)
        hold on
        
        % plot the histogram with a normal distribution fit
        h = histfit(data_rsEEG(:, m, p), 10, 'kernel');
        h(1).FaceColor = colours((m-1)*2 + 1, :); h(1).EdgeColor = [1 1 1];
        h(2).Color = [0 0 0]; h(2).LineWidth = 1;
        
        % get limits
        yl = get(gca, 'ylim'); xl = get(gca, 'xlim');
        
        % add mean of the distribution
        mu = mean(data_rsEEG(:, m, p));
        text(xl(1), yl(2)*-0.1, sprintf('mu = %1.3f', mu), 'FontSize', 10, 'FontWeight', 'bold')
        
        % add skewness and kurtosis
        s = skewness(data_rsEEG(:, m, p)); k = kurtosis(data_rsEEG(:, m, p));
        text(xl(1), yl(2)*-0.2, sprintf('%1.3f, %1.3f', s, k), 'FontSize', 10)
                       
        % set title
        if m == 1
            title(DV_rsEEG{p})
        end
    end
end   
clear p m h yl xl mu s k

% parameters
fig.Position = [500 400 750 500];
sgtitle('rsEEG - change in dependent variables')

% update counter
figure_counter = figure_counter + 1;

% QQ plot
fig = figure(figure_counter);
for p = 1:length(DV_rsEEG)
    for m = 1:length(medication)        
        % launch the plot
        subplot(length(medication), length(DV_rsEEG), (m-1)*length(DV_rsEEG) + p)
        hold on
        
        % plot the histogram with a normal distribution fit
        q = qqplot(data_rsEEG(:, m, p));
        title(''); xlabel(''); ylabel('')
        
        % get limits
        yl = get(gca, 'ylim'); xl = get(gca, 'xlim');
                       
        % set title
        if m == 1
            title(DV_rsEEG{p})
        end
    end
end   
clear p m q yl xl

% parameters
fig.Position = [200 400 1000 500];
sgtitle('rsEEG - change in dependent variables')

% update counter
figure_counter = figure_counter + 1;

%% 7) EFFECT OF ALPRAZOLAM - RS-EEG: export for R
GABA_YC_medication_rsEEG = table;
row_counter = 1;
for p = 1:length(participant) 
    for m = 1:length(medication)  
        % independent variables
        GABA_YC_medication_rsEEG.subject(row_counter) = participant(p);
        GABA_YC_medication_rsEEG.medication(row_counter) = medication(m);
                
        % dependent variables
        for k = 1:length(DV_rsEEG)             
            statement = ['GABA_YC_medication_rsEEG.' DV_rsEEG{k} '(row_counter) = data_rsEEG(p, m, k);'];
            eval(statement)             
        end
                    
        % update the counter
        row_counter = row_counter + 1;
    end
end
clear m p k row_counter
writetable(GABA_YC_medication_rsEEG, [folder_output '\GABA_YC_medication_rsEEG.csv'])

