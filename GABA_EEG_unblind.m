%% RS-EEG PREPROCESSING 2
% Written by Dominika for GABA-AD project (2021)
% 
% 1) Reads 'treatments.xlsx' and creates a cell array with medication
%    information for each dataset
% 2) Replaces session information with medication information, saves
% 3) Merges all subjects in group dataset 
%    for each category --> calls a lw6 function 'RLW_merge_epochs' 

%% parameters
clear all
clc

% dataset
% ----------------------------
target = 'AG';
prefix = 'avg bl icfilt ica26 crop but fft-notchfilt prefilt prea28 visual reref ds art-sup bl ep dc chan-select';
% ----------------------------
group = 'YC';
participant = 1:20;  
medication = {'placebo' 'alprazolam'};  
session = {'S1' 'S2'};      
time = {'pre' 'post'};
stimulus = {'CS 80' 'TS 120' 'ppTMS 80120'};
prefix_new = 'GABA';
stimulus_new = {'CS' 'TS' 'ppTMS'};

% load and adjust the medication info
med_session = table2cell(readtable('treatments.xlsx'));
for a = 1 : size(med_session, 1)
    for b = 1 : length(session)
        if strcmp(med_session{a, b+1}, 'zyrtec')
            med_session{a, b+1} = 'placebo'; 
        else
            med_session{a, b+1} = 'alprazolam'; 
        end
    end
end
clear a b

%% M1
% loop through datasets  
for p = 1:length(participant)
    % define participant
    if p < 10
        subj = ['0' num2str(participant(p))];
    else
        subj = num2str(participant(p));
    end
    
    % loop through datasets
    for s = 1:length(session)
        this_session = med_session(participant(p), s+1);
        for t = 1:length(time) 
            for c = 1:length(stimulus) 
                name_old = [prefix ' ' stimulus{c} ' ' group num2str(participant(p)) ' ' session{s} ' '  time{t} ' ' target];
                load([name_old '.mat']);
                load([name_old '.lw6'], '-mat');

                % save under new name
                name_new = [prefix_new ' ' group ' ' subj ' ' target ' ' char(this_session) ' '  time{t} ' ' stimulus_new{c}];
                header.name = name_new;
                save([name_new '.mat'], 'data');
                save([name_new '.lw6'], 'header');     
            end
        end
    end    
    message = ['Participant n.' num2str(participant(p)) ' finished.'];
    disp(message)
end
clear p s t c f name_old name_new data header

%% AG
% loop through datasets  
for p = 1:length(participant)
    % define participant
    if p < 10
        subj = ['0' num2str(participant(p))];
    else
        subj = num2str(participant(p));
    end
    
    % loop through datasets
    for s = 1:length(session)
        this_session = med_session(participant(p), s+1);
        for t = 1:length(time) 
            name_old = [prefix ' ' group num2str(participant(p)) ' ' session{s} ' ' target ' '  time{t}];
            load([name_old '.mat']);
            load([name_old '.lw6'], '-mat');

            % save under new name
            name_new = [prefix_new ' ' group ' ' subj ' ' target ' ' char(this_session) ' '  time{t}];
            header.name = name_new;
            save([name_new '.mat'], 'data');
            save([name_new '.lw6'], 'header');     
        end
    end    
    message = ['Participant n.' num2str(participant(p)) ' finished.'];
    disp(message)
end
clear p s t c f name_old name_new data header

%% merge epochs of each category
% create 'merge_idx' variable for epoch merging
merge_idx = 1:numel(participant);

% loop through conditions
for m = 1:numel(medication)
    for t = 1:numel(time)   
        for s = 1:length(stimulus)
            disp(['>>> ' medication{m} ': ' time{t} ' - ' stimulus_new{s} ' <<<'])    

            % create 'datasets' variable for epoch merging
            datasets = struct;
            for p = 1:numel(participant)
                % define participant
                if p < 10
                    subj = ['0' num2str(participant(p))];
                else
                    subj = num2str(participant(p));
                end

                % load deader and data
                load([prefix_new ' ' group ' ' subj ' ' target ' ' medication{m} ' '  time{t} ' ' stimulus_new{s} '.lw6'], '-mat');
                datasets(p).header = header;
                load([prefix_new ' ' group ' ' subj ' ' target ' ' medication{m} ' '  time{t} ' ' stimulus_new{s} '.mat']);
                datasets(p).data = data;
            end
            clear data header

            % merge datasets
            [header,data,message_string] = RLW_merge_epochs(datasets,merge_idx);                

            % save the data
            name_new = ['merged ' prefix_new ' ' group ' ' target ' ' medication{m} ' '  time{t} ' ' stimulus_new{s}];
            save([name_new '.mat'], 'data')

            % save the header
            header.name = name_new;
            save([name_new '.lw6'], 'header')      
        end
    end    
end
clear data header name_new name_old datasets message_string
disp('>>> All done. <<<')

