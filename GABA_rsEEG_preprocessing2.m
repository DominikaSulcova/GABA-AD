%% RS-EEG PREPROCESSING 2
% Written by Dominika for GABA-AD project (2021)
% 
% 1) Reads 'treatments.xlsx' and creates a cell array with medication
%    information for each dataset
% 2) Replaces session information with medication information, saves
% 3) After averaging (done in letswave), merges all subjects in group dataset 
%    for each category --> calls a lw6 function 'RLW_merge_epochs' 

%% parameters
clear all
clc

% dataset
participant = 1:20;  
medication = {'placebo' 'alprazolam'};  
session = {'S1' 'S2'};      
time = {'pre' 'post'};
condition = {'open' 'closed'};
frequency = {'low' 'high'};
prefix_1 = 'han icfilt assigned ica'; 
prefix_2 = 'crop_120 45Hz_but4 but reref chan-select ds YC';

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
            for c = 1:length(condition)
                for f = 1:length(frequency)
                    % rename fully filtered data (all artifactual ICAs removed)
                    name_old = ['fft ' prefix_1 ' ' frequency{f} ' ' prefix_2 num2str(participant(p)) ' ' session{s} ' ' time{t} ' EEGcont ' condition{c}];
                    load([name_old '.mat']);
                    load([name_old '.lw6'], '-mat');
                    
                    % save under new name
                    name_new = ['fft ' prefix_1 ' ' frequency{f} ' ' prefix_2 num2str(participant(p)) ' ' char(this_session) ' ' time{t} ' EEGcont ' condition{c}];
                    header.name = name_new;
                    save([name_new '.mat'], 'data');
                    save([name_new '.lw6'], 'header');                        
                end 
            end
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
        for c = 1:numel(condition)
            disp(['>>> ' medication{m} ' - ' time{t} ' - ' condition{c} ' <<<'])
            for f = 1:numel(frequency)
                disp(['..... ' frequency{f} ' frequency .....'])              
                
                % create 'datasets' variable for epoch merging
                datasets = struct;
                for p = 1:numel(participant)
                    load(['avg fft ' prefix_1 ' ' frequency{f} ' ' prefix_2 num2str(participant(p)) ' ' medication{m} ' ' time{t} ' EEGcont ' condition{c} '.lw6'], '-mat');
                    datasets(p).header = header;
                    load(['avg fft ' prefix_1 ' ' frequency{f} ' ' prefix_2 num2str(participant(p)) ' ' medication{m} ' ' time{t} ' EEGcont ' condition{c} '.mat']);
                    datasets(p).data = data;
                end
                clear data header

                % merge datasets
                [header,data,message_string] = RLW_merge_epochs(datasets,merge_idx);                
                
                % save the data
                name_new = ['merged fft ' frequency{f} ' ' medication{m} ' ' time{t} ' EEGcont ' condition{c}];
                save([name_new '.mat'], 'data')

                % save the header
                header.name = name_new;
                save([name_new '.lw6'], 'header')                        
            end
        end
    end    
end
clear data header name_new name_old datasets message_string
disp('>>> All done. <<<')

