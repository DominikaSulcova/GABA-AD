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
target = 'AG';
participant = 1:20;  
medication = {'placebo' 'alprazolam'};  
session = {'S1' 'S2'};      
time = {'pre' 'post'};
prefix = 'avg bl icfilt ica26 crop but fft-notchfilt prefilt prea28 visual reref ds art-sup bl ep dc chan-select YC';

% load and adjust the medication info
med_session = table2cell(readtable('treatments.xlsx'));
for a = 1 : size(med_session, 1)
    for b = 1 : length(session)
        if med_session{a, b+1}(1) == 'z'
            med_session{a, b+1} = 'placebo'; 
        else
            med_session{a, b+1} = 'alprazolam'; 
        end
    end
end
clear a b

%% replace the session number in the name by corresponding medication
% loop through datasets  
for p = 1:length(participant)
    for s = 1:length(session)
        this_session = med_session(participant(p), s+1);
        for t = 1:length(time)       
            name_old = [prefix num2str(participant(p)) ' ' session{s} ' ' target ' '  time{t}];
            load([name_old '.mat']);
            load([name_old '.lw6'], '-mat');

            % save under new name
            name_new = [prefix num2str(participant(p)) ' ' char(this_session) ' ' target ' '  time{t}];
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
        disp(['>>> ' medication{m} ' - ' time{t}  ' <<<'])          
        % create 'datasets' variable for epoch merging
        datasets = struct;
        for p = 1:numel(participant)
            load([prefix num2str(participant(p)) ' ' medication{m} ' ' target ' '  time{t} '.lw6'], '-mat');
            datasets(p).header = header;
            load([prefix num2str(participant(p)) ' ' medication{m} ' ' target ' '  time{t} '.mat']);
            datasets(p).data = data;
        end
        clear data header

        % merge datasets
        [header,data,message_string] = RLW_merge_epochs(datasets,merge_idx);                

        % save the data
        name_new = ['merged ' medication{m} ' ' time{t} ' AG'];
        save([name_new '.mat'], 'data')

        % save the header
        header.name = name_new;
        save([name_new '.lw6'], 'header')                        
    end    
end
clear data header name_new name_old datasets message_string
disp('>>> All done. <<<')

